// cabart_ess.cpp -- AFT-BART (censored log-normal) that emits per-leaf
// summaries during the in-stream Gibbs sweep, so the R driver can
// compute leaf-level ELIR ESS post-sweep via RBesT_lite.  Survival
// counterpart of cwbart_ess.cpp.
//
// Algorithm.  At each iteration b, for each tree h in turn:
//   1. Compute partial residual r[k] = z[k] - sum_{h' != h} g_{h'}(x_k),
//      where z[k] is the *latent log-time*: observed iy[k] = log(Y_obs)
//      for events (delta==1), or last iteration's rtnorm-imputed value
//      for censored observations (delta==0).
//   2. Update tree h via heterbd + heterdrmu (heteroscedastic because
//      svec[k] = iw[k] * sigma carries per-obs weights).
//   3. EMIT calibration hook: walk leaves of t[h], compute ybar and n
//      per leaf, append a row (b, h, leaf_id, ybar, n, sigma_prev) to
//      the output table.  sigma_prev is the sigma in scope during this
//      sweep (sigma^(b-1) in the algorithm).
// At the end of the per-tree sweep:
//   - Draw a new sigma from its conditional InvGamma posterior
//     (rss-based, using all uncensored + currently-imputed latent
//     residuals just as cabart.cpp does).
//   - For censored k: re-impute z[k] ~ TruncNorm(bm.f(k), iy[k], svec[k]).
//
// This matches cabart.cpp's main loop exactly, with two changes:
//   - draw() is split into per-tree pieces so the calibration hook can
//     fire between drmu() and the re-fit.
//   - We only support type=1 (continuous AFT outcome); pbart/lbart
//     branches from cabart.cpp are omitted as they are unused here.

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

typedef std::vector<double> v1d;

using std::cout;
using std::endl;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Header order matches aBART/cabart.cpp.
#include "include_ess/rn.h"
#include "include_ess/tree.h"
#include "include_ess/treefuns.h"
#include "include_ess/info.h"
#include "include_ess/bartfuns.h"
#include "include_ess/bart.h"
#include "include_ess/heterbartfuns.h"
#include "include_ess/heterbd.h"
#include "include_ess/heterbart.h"
#include "include_ess/rtnorm.h"

//-------------------------------------------------------------------
// calibrate_tree: walk the leaves of t[j], compute (ybar, n) from the
// partial residual r[k] (exposed via bm.getr), and append one row per
// leaf to `rows`.  Each row is six doubles:
//   (b, h, leaf_id, ybar, n, sigma_prev)
//-------------------------------------------------------------------
static void calibrate_tree(heterbart& bm, size_t b, size_t j, double sigma_prev,
                            v1d& rows)
{
   size_t n     = bm.getn_();
   xinfo& xi    = bm.get_xi();
   tree&  treej = bm.gettree(j);

   tree::npv leaves;
   treej.getbots(leaves);
   size_t L = leaves.size();
   if(L == 0) return;

   v1d leaf_sum(L, 0.0);
   std::vector<size_t> leaf_n(L, 0);

   for(size_t k=0;k<n;k++) {
      tree::tree_p bn = treej.bn(bm.getx_row(k), xi);
      size_t li = 0;
      for(; li<L; ++li) if(leaves[li] == bn) break;
      if(li == L) continue;
      leaf_sum[li] += bm.getr(k);
      leaf_n[li]   += 1;
   }

   for(size_t li=0; li<L; ++li) {
      size_t nl = leaf_n[li];
      if(nl == 0) continue;
      double ybar = leaf_sum[li] / (double)nl;
      rows.push_back((double)b);
      rows.push_back((double)j);
      rows.push_back((double)li);
      rows.push_back(ybar);
      rows.push_back((double)nl);
      rows.push_back(sigma_prev);
   }
}

//===================================================================
// Rcpp entry point.  Mirrors cwbart_ess (Gaussian) but uses heterbart
// + per-iteration rtnorm imputation of censored latent log-times.
//
// Input convention (matches cabart.cpp):
//   iy[k]   = log(Y_obs_k) for both event and censored observations
//   delta[k] = 1 if event, 0 if censored
//   iw[k]   = per-obs weight (usually all 1.0; passed through to svec).
//===================================================================
// [[Rcpp::export]]
Rcpp::List cabart_ess(
    int                 in_,
    int                 ip,
    Rcpp::NumericVector ix,           // x: column-stacked p x n
    Rcpp::NumericVector iy,           // log(Y_obs)
    Rcpp::IntegerVector idelta,       // 1 = event, 0 = censored
    int                 im,
    Rcpp::IntegerVector inc,
    int                 ind,
    int                 iburn,
    double              ipower,
    double              ibase,
    double              itau,
    double              inu,
    double              ilambda,
    double              isigest,
    Rcpp::NumericVector iiw,          // per-obs weights w[k]; pass rep(1, n)
    int                 iverbose,
    Rcpp::Nullable<Rcpp::NumericMatrix> iXinfo = R_NilValue
)
{
   using namespace Rcpp;

   size_t n     = (size_t)in_;
   size_t p     = (size_t)ip;
   size_t m     = (size_t)im;
   size_t nd    = (size_t)ind;
   size_t burn  = (size_t)iburn;
   double mybeta = ipower;
   double alpha  = ibase;
   double tau    = itau;
   double nu     = inu;
   double lambda = ilambda;
   double sigma  = isigest;
   bool   verbose = (iverbose != 0);

   double* ix_p   = REAL(ix);
   double* iy_p   = REAL(iy);
   int*    delta  = INTEGER(idelta);
   double* iw     = REAL(iiw);
   int*    nc     = INTEGER(inc);

   //------ Latent z[k] and per-obs SD svec[k].
   // z[k] starts at iy[k] (the observed log-time); censored rows get
   // re-imputed via rtnorm at the end of each sweep.
   std::vector<double> z(n), svec(n);
   for(size_t k=0;k<n;k++) {
      z[k]    = iy_p[k];
      svec[k] = iw[k] * sigma;
   }

   //------ AFT-BART (heterbart) setup, mirroring cabart.cpp.
   heterbart bm(m);
   bm.setprior(alpha, mybeta, tau);
   bm.setdata(p, n, ix_p, &z[0], nc);
   if(iXinfo.isNotNull()) {
      NumericMatrix Xi(iXinfo.get());
      xinfo xiv;
      xiv.resize(p);
      for(size_t i=0;i<p;++i) {
         xiv[i].resize(nc[i]);
         for(int j=0;j<nc[i];++j) xiv[i][j] = Xi(i, j);
      }
      bm.setxinfo(xiv);
   }

   GetRNGstate();
   arn gen;

   //------ Output buffers.
   v1d rows;
   rows.reserve(nd * m * 8 * 6);
   NumericVector sigma_draws(nd);

   //------ Gibbs sweeps.
   double df = (double)n + nu;
   double sigma_prev = sigma;
   for(size_t i=0;i<(nd+burn);++i) {
      // Per-tree sweep with calibration emission after heterdrmu().
      for(size_t j=0;j<m;++j) {
         bm.draw_per_tree(j, &svec[0], gen);
         if(i >= burn) {
            calibrate_tree(bm, i - burn, j, sigma_prev, rows);
         }
         bm.finish_per_tree(j);
      }

      // AFT-BART residual-scale update.  rss uses the *current* latent
      // log-times z[k] (which include the rtnorm imputations of the
      // previous sweep for censored rows).
      double rss = 0.0;
      for(size_t k=0;k<n;++k) {
         double e = (z[k] - bm.f(k)) / iw[k];
         rss += e * e;
      }
      sigma_prev = std::sqrt((nu * lambda + rss) / gen.chi_square(df));

      // Re-impute censored latent log-times: z[k] ~ N(bm.f(k), svec[k])
      // truncated below by iy[k] = log(Y_obs).
      for(size_t k=0;k<n;++k) {
         svec[k] = iw[k] * sigma_prev;
         if(delta[k] == 0) z[k] = rtnorm(bm.f(k), iy_p[k], svec[k], gen);
      }

      if(i >= burn) sigma_draws[i - burn] = sigma_prev;

      if(verbose && ((i + 1) % 100 == 0)) {
         Rcpp::Rcout << "  iter " << (i+1) << "/" << (nd+burn)
                      << (i < burn ? " (warmup)" : "")
                      << "  sigma=" << sigma_prev << endl;
      }
   }

   PutRNGstate();

   //------ Pack output.
   size_t nrows = rows.size() / 6;
   NumericMatrix leaf_table((int)nrows, 6);
   for(size_t r=0; r<nrows; ++r)
      for(size_t c=0; c<6; ++c)
         leaf_table(r, c) = rows[r*6 + c];
   colnames(leaf_table) = CharacterVector::create("b","h","leaf","ybar","n","sigma");

   List ret;
   ret["leaf_table"]  = leaf_table;
   ret["sigma_draws"] = sigma_draws;
   ret["nd"]          = (int)nd;
   ret["burn"]        = (int)burn;
   return ret;
}
