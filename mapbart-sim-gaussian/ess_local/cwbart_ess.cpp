// cwbart_ess.cpp -- Gaussian-outcome standard BART that emits per-leaf
// partial-residual summaries during the in-stream Gibbs sweep, so the
// R driver can compute leaf-level ELIR ESS post-sweep via RBesT_lite.
//
// Implementation note.  Algorithm 1 of main.tex specifies the BART
// partial residual r_i^(h,b) for tree h at iteration b, formed from
// iteration-b leaf values for h' < h and iteration-(b-1) values for
// h' > h.  This is exactly what BART's internal Gibbs sweep already
// computes (see wBART/include.w/bart.h, void bart::draw -- the array
// r[k] is the partial residual at the moment tree j is being updated).
// To expose r[k] between drmu() and the re-fit, the forked header
// include_ess/bart.h splits draw() into draw_per_tree(j) +
// finish_per_tree(j), and exposes accessors getr/getx_row/get_xi etc.
// This file's calibrate_tree() runs in between and emits one row per
// leaf to the output table.
//
// No ELIR / gMAP / mixture math happens in this file -- those run in R.

#include <Rcpp.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
typedef std::vector<double> v1d;

using std::cout;
using std::endl;

// Header order matches wBART/cwbart.cpp.  Default RNG path (YesRcpp)
// uses R::norm_rand etc and lives inside Rcpp's RNGScope, so seeding
// is via set.seed() on the R side before this function is called.
#include "include_ess/rn.h"
#include "include_ess/tree.h"
#include "include_ess/treefuns.h"
#include "include_ess/info.h"
#include "include_ess/bartfuns.h"
#include "include_ess/bd.h"
#include "include_ess/bart.h"

//-------------------------------------------------------------------
// calibrate_tree: walk the leaves of t[j], compute (ybar, n) from the
// BART partial residual r[k] (exposed via bm.getr), and append one row
// per leaf to `rows`.  Each row is six doubles:
//   (b, h, leaf_id, ybar, n, sigma_prev)
// where sigma_prev is the sigma in scope during this sweep, captured
// BEFORE the post-tree-sweep sigma update.
//-------------------------------------------------------------------
static void calibrate_tree(bart& bm, size_t b, size_t j, double sigma_prev,
                            v1d& rows)
{
   size_t n     = bm.getn_();
   xinfo& xi    = bm.get_xi();
   tree&  treej = bm.gettree(j);

   tree::npv leaves;
   treej.getbots(leaves);
   size_t L = leaves.size();
   if(L == 0) return;

   // Per-leaf accumulators (small: L is typically <= 20 under BART
   // tree prior at this scale).
   v1d leaf_sum(L, 0.0);
   std::vector<size_t> leaf_n(L, 0);

   // Walk each observation, find which leaf of t[j] it lands in,
   // accumulate r[k].
   for(size_t k=0;k<n;k++) {
      tree::tree_p bn = treej.bn(bm.getx_row(k), xi);
      size_t li = 0;
      for(; li<L; ++li) if(leaves[li] == bn) break;
      if(li == L) continue;   // shouldn't happen
      leaf_sum[li] += bm.getr(k);
      leaf_n[li]   += 1;
   }

   for(size_t li=0; li<L; ++li) {
      size_t nl = leaf_n[li];
      if(nl == 0) continue;   // empty leaf -- not expected under n_min>=5
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
// Rcpp entry point: returns a NumericMatrix with columns
//   [b, h, leaf_id, ybar, n, sigma]
// and as many rows as total leaves across all (b, h) cells.
//===================================================================
// [[Rcpp::export]]
Rcpp::List cwbart_ess(
    int                 in_,
    int                 ip,
    Rcpp::NumericVector ix,           // x: column-stacked p x n (pre-transposed in R)
    Rcpp::NumericVector iy,           // y: length-n training response
    int                 im,           // number of trees H
    Rcpp::IntegerVector inc,          // numcut per predictor (length p)
    int                 ind,          // B: number of post-warmup retained draws
    int                 iburn,        // warmup
    double              ipower,       // BART tree-depth beta
    double              ibase,        // BART tree-depth alpha
    double              itau,         // BART leaf-prior tau
    double              inu,          // BART sigma-prior nu
    double              ilambda,      // BART sigma-prior lambda
    double              isigest,      // initial sigma
    int                 iverbose,     // 1 = print progress every 100 iters
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

   double* ix_p = REAL(ix);
   double* iy_p = REAL(iy);
   int*    nc   = INTEGER(inc);

   //------ BART setup (mirrors cwbart.cpp) ------
   bart bm(m);
   bm.setprior(alpha, mybeta, tau);
   bm.setdata(p, n, ix_p, iy_p, nc);
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

   //------ RNG (uses R's RNGScope; set.seed() on R side controls it) ------
   GetRNGstate();
   arn gen;

   //------ Output buffers ------
   v1d rows;
   rows.reserve(nd * m * 8 * 6);   // rough pre-allocation
   NumericVector sigma_draws(nd);
   NumericVector accept_draws(nd);

   //------ Gibbs sweeps ------
   double rss, restemp;
   double sigma_prev = sigma;
   for(size_t i=0;i<(nd+burn);++i) {
      // Per-tree sweep with calibration emission after drmu()
      size_t accept_count = 0;
      for(size_t j=0;j<m;++j) {
         bm.draw_per_tree(j, sigma_prev, gen);
         if(i >= burn) {
            calibrate_tree(bm, i - burn, j, sigma_prev, rows);
         }
         bm.finish_per_tree(j);
      }
      // BART residual-scale update (post-tree-sweep).  sigma_prev for
      // the *next* sweep becomes the just-sampled sigma.
      rss = 0.0;
      for(size_t k=0;k<n;++k) {
         restemp = iy_p[k] - bm.f(k);
         rss += restemp * restemp;
      }
      sigma_prev = std::sqrt((nu * lambda + rss) / gen.chi_square(n + nu));

      if(i >= burn) sigma_draws[i - burn] = sigma_prev;

      if(verbose && ((i + 1) % 100 == 0)) {
         Rcpp::Rcout << "  iter " << (i+1) << "/" << (nd+burn)
                      << (i < burn ? " (warmup)" : "")
                      << "  sigma=" << sigma_prev << endl;
      }
   }

   PutRNGstate();

   //------ Pack output ------
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
