// lrcbart: locally robust commensurate BART, Rcpp core.
//
// Model (method_spec.md v2, Section A1):
//   y_i = f(x_i) + D_i g(x_i) + eps_i,  eps_i ~ N(0, sigma_{src_i}^2),
//   D_i = 1{RCT control}, src = 1 (RCT control) or 2 (RWD control).
// f: standard BART ensemble on the pooled controls, leaf prior N(0, lambda_f^2).
// g: sparse ensemble on the RCT controls, spike-and-slab leaf prior
//   theta | z ~ z N(0, tau0^2) + (1 - z) N(0, tau1^2), z ~ Bern(w).
//
// One Tree class is shared by both ensembles; the tree-move code is written
// once and dispatches to the f or g leaf marginal through the Ensemble flag.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>

using namespace Rcpp;

namespace {

const double LOG2PI = 1.8378770664093453;
const double NEG_INF = -std::numeric_limits<double>::infinity();

typedef std::vector<std::vector<double> > CutPoints;

inline double ldnorm(double x, double m, double v) {
  double d = x - m;
  return -0.5 * (LOG2PI + std::log(v) + d * d / v);
}

inline double lse2(double a, double b) {
  if (a == NEG_INF) return b;
  if (b == NEG_INF) return a;
  double m = a > b ? a : b;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// ---------------------------------------------------------------------------
// Leaf marginals (spec Section B(i)), all in log space.
// ---------------------------------------------------------------------------

// f leaf: heteroscedastic Chipman factor.
// A = n1/sig1^2 + n2/sig2^2, B = S1/sig1^2 + S2/sig2^2.
inline double f_leaf_lm(double n1, double n2, double s1, double s2,
                        double sig1sq, double sig2sq, double lam2) {
  double A = n1 / sig1sq + n2 / sig2sq;
  double B = s1 / sig1sq + s2 / sig2sq;
  double c = 1.0 + lam2 * A;
  return -0.5 * std::log(c) + lam2 * B * B / (2.0 * c);
}

// g leaf: two-component mixture in ybar1, divided by the reference N(ybar; 0, v1).
inline double g_leaf_lm(double n1, double s1, double sig1sq,
                        double tau0sq, double tau1sq, double w) {
  if (n1 <= 0) return 0.0;
  double yb = s1 / n1, v1 = sig1sq / n1;
  double l0 = (w > 0.0) ? std::log(w) + ldnorm(yb, 0.0, tau0sq + v1) : NEG_INF;
  double l1 = (w < 1.0) ? std::log(1.0 - w) + ldnorm(yb, 0.0, tau1sq + v1) : NEG_INF;
  return lse2(l0, l1) - ldnorm(yb, 0.0, v1);
}

// P(z = 1 | ybar1) (spec B(iii)); n1 = 0 gives the prior w.
inline double g_leaf_pspike(double n1, double s1, double sig1sq,
                            double tau0sq, double tau1sq, double w) {
  if (w >= 1.0) return 1.0;
  if (w <= 0.0) return 0.0;
  if (n1 <= 0) return w;
  double yb = s1 / n1, v1 = sig1sq / n1;
  double l0 = std::log(w) + ldnorm(yb, 0.0, tau0sq + v1);
  double l1 = std::log(1.0 - w) + ldnorm(yb, 0.0, tau1sq + v1);
  return std::exp(l0 - lse2(l0, l1));
}

// Collapsed block draw of (z, theta) for a g leaf.
inline void g_leaf_draw(double n1, double s1, double sig1sq,
                        double tau0sq, double tau1sq, double w,
                        int& z, double& theta) {
  double p1 = g_leaf_pspike(n1, s1, sig1sq, tau0sq, tau1sq, w);
  z = (unif_rand() < p1) ? 1 : 0;
  double tz = z ? tau0sq : tau1sq;
  if (n1 <= 0) {
    theta = std::sqrt(tz) * norm_rand();
    return;
  }
  double yb = s1 / n1, v1 = sig1sq / n1;
  double m = tz * yb / (tz + v1);
  double V = tz * v1 / (tz + v1);
  theta = m + std::sqrt(V) * norm_rand();
}

// f leaf draw: N(lam2 B / (1 + lam2 A), lam2 / (1 + lam2 A)).
inline double f_leaf_draw(double n1, double n2, double s1, double s2,
                          double sig1sq, double sig2sq, double lam2) {
  double A = n1 / sig1sq + n2 / sig2sq;
  double B = s1 / sig1sq + s2 / sig2sq;
  double c = 1.0 + lam2 * A;
  return lam2 * B / c + std::sqrt(lam2 / c) * norm_rand();
}

// Truncated normal N(mu, sigma^2) restricted to (a, Inf). Robert (1995)
// exponential rejection in the tail, plain rejection otherwise.
double rtnorm_lower(double mu, double sigma, double a) {
  double as = (a - mu) / sigma;
  double z;
  if (as < 0.5) {
    do { z = norm_rand(); } while (z < as);
  } else {
    double lam = 0.5 * (as + std::sqrt(as * as + 4.0));
    for (;;) {
      z = as + exp_rand() / lam;
      double rho = std::exp(-0.5 * (z - lam) * (z - lam));
      if (unif_rand() < rho) break;
    }
  }
  return mu + sigma * z;
}

// scaled-Inv-chi^2(nu, s^2) draw: nu s^2 / chi^2_nu.
inline double rsinvchisq(double nu, double ssq) {
  return nu * ssq / R::rchisq(nu);
}

// scaled-Inv-chi^2(nu, s^2) truncated to (0, upper): chi^2 > nu s^2 / upper.
inline double rsinvchisq_upper(double nu, double ssq, double upper) {
  double c = nu * ssq / upper;
  double pu = R::pchisq(c, nu, 0, 0);  // P(chi2 > c)
  if (pu <= 0.0) return upper * (1.0 - 1e-12);  // numerically at the bound
  double u = unif_rand() * pu;
  double q = R::qchisq(u, nu, 0, 0);
  return nu * ssq / q;
}

// scaled-Inv-chi^2(nu, s^2) truncated to (lower, Inf): chi^2 < nu s^2 / lower.
inline double rsinvchisq_lower(double nu, double ssq, double lower) {
  double c = nu * ssq / lower;
  double pl = R::pchisq(c, nu, 1, 0);  // P(chi2 < c)
  if (pl <= 0.0) return lower * (1.0 + 1e-12);
  double u = unif_rand() * pl;
  double q = R::qchisq(u, nu, 1, 0);
  return nu * ssq / q;
}

// ---------------------------------------------------------------------------
// Tree
// ---------------------------------------------------------------------------

struct Node {
  int var, cut, left, right, parent;
  double mu;
  int z;
};

struct Tree {
  std::vector<Node> nd;
  Tree() {
    Node r; r.var = -1; r.cut = -1; r.left = -1; r.right = -1; r.parent = -1;
    r.mu = 0.0; r.z = 1;
    nd.push_back(r);
  }
  bool leaf(int i) const { return nd[i].left < 0; }
  bool stump() const { return nd.size() == 1; }
  int depth(int i) const {
    int d = 0;
    while (nd[i].parent >= 0) { i = nd[i].parent; ++d; }
    return d;
  }
  void leaves(std::vector<int>& out) const {
    out.clear();
    for (size_t i = 0; i < nd.size(); ++i) if (leaf((int)i)) out.push_back((int)i);
  }
  // internal nodes whose two children are both leaves
  void nogs(std::vector<int>& out) const {
    out.clear();
    for (size_t i = 0; i < nd.size(); ++i)
      if (!leaf((int)i) && leaf(nd[i].left) && leaf(nd[i].right)) out.push_back((int)i);
  }
  int nleaves() const { int k = 0; for (size_t i = 0; i < nd.size(); ++i) if (leaf((int)i)) ++k; return k; }
  int nnogs() const { std::vector<int> v; nogs(v); return (int)v.size(); }
  // cut indices for variable v still available at node i given its ancestors
  void bounds(int i, int v, int ncut, int& lo, int& hi) const {
    lo = 0; hi = ncut - 1;
    int c = i;
    while (nd[c].parent >= 0) {
      int p = nd[c].parent;
      if (nd[p].var == v) {
        if (nd[p].left == c) hi = std::min(hi, nd[p].cut - 1);
        else lo = std::max(lo, nd[p].cut + 1);
      }
      c = p;
    }
  }
  void grow(int l, int v, int c) {
    Node a; a.var = -1; a.cut = -1; a.left = -1; a.right = -1; a.parent = l; a.mu = 0.0; a.z = 1;
    nd[l].var = v; nd[l].cut = c;
    nd[l].left = (int)nd.size(); nd.push_back(a);
    nd[l].right = (int)nd.size(); nd.push_back(a);
  }
  Tree pruned(int nog) const {
    std::vector<int> map(nd.size(), -1);
    Tree t; t.nd.clear();
    int l = nd[nog].left, r = nd[nog].right;
    for (size_t i = 0; i < nd.size(); ++i) {
      if ((int)i == l || (int)i == r) continue;
      map[i] = (int)t.nd.size();
      t.nd.push_back(nd[i]);
    }
    for (size_t i = 0; i < t.nd.size(); ++i) {
      Node& q = t.nd[i];
      if (q.parent >= 0) q.parent = map[q.parent];
      q.left = (q.left >= 0) ? map[q.left] : -1;
      q.right = (q.right >= 0) ? map[q.right] : -1;
      if (q.left < 0) { q.var = -1; q.cut = -1; }
    }
    return t;
  }
};

inline int find_leaf(const Tree& t, const std::vector<double>& Xv, int n, int i,
                     const CutPoints& cp) {
  int k = 0;
  while (t.nd[k].left >= 0) {
    const Node& q = t.nd[k];
    k = (Xv[(size_t)q.var * n + i] < cp[q.var][q.cut]) ? q.left : q.right;
  }
  return k;
}

// available split variables at node i
inline void avail_vars(const Tree& t, int i, const CutPoints& cp, std::vector<int>& out) {
  out.clear();
  for (size_t v = 0; v < cp.size(); ++v) {
    if (cp[v].empty()) continue;
    int lo, hi;
    t.bounds(i, (int)v, (int)cp[v].size(), lo, hi);
    if (lo <= hi) out.push_back((int)v);
  }
}

// ---------------------------------------------------------------------------
// Ensemble and sampler
// ---------------------------------------------------------------------------

enum Move { GROW = 0, PRUNE = 1, CHANGE = 2 };

struct Ensemble {
  bool isg;
  int H;
  double alpha, beta;
  int nmin;
  std::vector<Tree> trees;
  std::vector<std::vector<double> > fit;  // H x n
  std::vector<double> tot;               // n
  long acc[3], prop[3];
};

class Sampler {
public:
  int n, p;
  std::vector<double> Xv;   // column-major n x p
  std::vector<double> y;
  std::vector<int> src;     // 1 or 2
  std::vector<int> status;  // 1 observed, 0 censored (AFT); all 1 for Gaussian
  std::vector<double> logC; // censoring bound (log scale), used if status == 0
  CutPoints cp;
  int n1obs, n2obs;
  bool single_arm;

  // parameters
  double sig1sq, sig2sq, tau0sq, tau1sq, w, lam2;
  // hyperparameters
  double nu_sig, lam_sig1, lam_sig2, nu0, s0sq, nu1, s1sq, aw, bw;
  bool w_fixed, upd_tau0, upd_tau1, aft;
  double pg, pp, pc;

  Ensemble F, G;

  // buffers
  std::vector<int> leaf_of;
  std::vector<double> bn1, bn2, bs1, bs2;
  std::vector<double> resid;
  std::vector<int> ibuf, ibuf2;

  double pmove(const Tree& t, int m) const {
    if (t.stump()) return (m == GROW) ? 1.0 : 0.0;
    return (m == GROW) ? pg : (m == PRUNE ? pp : pc);
  }

  int pick_move(const Tree& t) const {
    if (t.stump()) return GROW;
    double u = unif_rand();
    if (u < pg) return GROW;
    if (u < pg + pp) return PRUNE;
    return CHANGE;
  }

  // assign all observations to leaves of t and accumulate leaf statistics
  void stats(const Tree& t, const Ensemble& E) {
    size_t m = t.nd.size();
    bn1.assign(m, 0.0); bn2.assign(m, 0.0); bs1.assign(m, 0.0); bs2.assign(m, 0.0);
    for (int i = 0; i < n; ++i) {
      int l = find_leaf(t, Xv, n, i, cp);
      leaf_of[i] = l;
      if (src[i] == 1) { bn1[l] += 1.0; bs1[l] += resid[i]; }
      else if (!E.isg) { bn2[l] += 1.0; bs2[l] += resid[i]; }
    }
  }

  double tree_loglik(const Tree& t, const Ensemble& E, bool check_min) {
    stats(t, E);
    double ll = 0.0;
    for (size_t l = 0; l < t.nd.size(); ++l) {
      if (!t.leaf((int)l)) continue;
      double cnt = E.isg ? bn1[l] : bn1[l] + bn2[l];
      if (check_min && cnt < E.nmin) return NEG_INF;
      if (E.isg) ll += g_leaf_lm(bn1[l], bs1[l], sig1sq, tau0sq, tau1sq, w);
      else ll += f_leaf_lm(bn1[l], bn2[l], bs1[l], bs2[l], sig1sq, sig2sq, lam2);
    }
    return ll;
  }

  double log_split_prior_ratio(const Ensemble& E, int d) const {
    // P(split at depth d) * P(no split at d+1)^2 / P(no split at d)
    double ps_d = E.alpha * std::pow(1.0 + d, -E.beta);
    double ps_d1 = E.alpha * std::pow(2.0 + d, -E.beta);
    return std::log(ps_d) + 2.0 * std::log1p(-ps_d1) - std::log1p(-ps_d);
  }

  void compute_resid(const Ensemble& E, int h) {
    if (E.isg) {
      for (int i = 0; i < n; ++i)
        resid[i] = (src[i] == 1) ? y[i] - F.tot[i] - G.tot[i] + E.fit[h][i] : 0.0;
    } else {
      for (int i = 0; i < n; ++i) {
        double gpart = (src[i] == 1 && G.H > 0) ? G.tot[i] : 0.0;
        resid[i] = y[i] - F.tot[i] + E.fit[h][i] - gpart;
      }
    }
  }

  void update_tree(Ensemble& E, int h) {
    compute_resid(E, h);
    Tree& t = E.trees[h];
    int mv = pick_move(t);
    E.prop[mv]++;
    if (mv == GROW) {
      t.leaves(ibuf);
      int l = ibuf[(int)(unif_rand() * ibuf.size())];
      avail_vars(t, l, cp, ibuf2);
      if (!ibuf2.empty()) {
        int v = ibuf2[(int)(unif_rand() * ibuf2.size())];
        int lo, hi;
        t.bounds(l, v, (int)cp[v].size(), lo, hi);
        int c = lo + (int)(unif_rand() * (hi - lo + 1));
        Tree tn = t;
        tn.grow(l, v, c);
        double lln = tree_loglik(tn, E, true);
        if (lln > NEG_INF) {
          double llc = tree_loglik(t, E, false);
          int d = t.depth(l);
          double lr = lln - llc + log_split_prior_ratio(E, d)
            + std::log(pmove(tn, PRUNE)) - std::log(pmove(t, GROW))
            + std::log((double)t.nleaves()) - std::log((double)tn.nnogs());
          if (std::log(unif_rand()) < lr) { t = tn; E.acc[GROW]++; }
        }
      }
    } else if (mv == PRUNE) {
      t.nogs(ibuf);
      int nog = ibuf[(int)(unif_rand() * ibuf.size())];
      Tree tn = t.pruned(nog);
      double lln = tree_loglik(tn, E, false);
      double llc = tree_loglik(t, E, false);
      int d = t.depth(nog);
      double lr = lln - llc - log_split_prior_ratio(E, d)
        + std::log(pmove(tn, GROW)) - std::log(pmove(t, PRUNE))
        + std::log((double)t.nnogs()) - std::log((double)tn.nleaves());
      if (std::log(unif_rand()) < lr) { t = tn; E.acc[PRUNE]++; }
    } else {
      t.nogs(ibuf);
      int nog = ibuf[(int)(unif_rand() * ibuf.size())];
      avail_vars(t, nog, cp, ibuf2);
      if (!ibuf2.empty()) {
        int v = ibuf2[(int)(unif_rand() * ibuf2.size())];
        int lo, hi;
        t.bounds(nog, v, (int)cp[v].size(), lo, hi);
        int c = lo + (int)(unif_rand() * (hi - lo + 1));
        Tree tn = t;
        tn.nd[nog].var = v; tn.nd[nog].cut = c;
        double lln = tree_loglik(tn, E, true);
        if (lln > NEG_INF) {
          double llc = tree_loglik(t, E, false);
          double lr = lln - llc;
          if (std::log(unif_rand()) < lr) { t = tn; E.acc[CHANGE]++; }
        }
      }
    }
    // draw leaf parameters for the (possibly updated) tree
    stats(t, E);
    for (size_t l = 0; l < t.nd.size(); ++l) {
      if (!t.leaf((int)l)) continue;
      if (E.isg) {
        int z; double th;
        g_leaf_draw(bn1[l], bs1[l], sig1sq, tau0sq, tau1sq, w, z, th);
        t.nd[l].z = z; t.nd[l].mu = th;
      } else {
        t.nd[l].mu = f_leaf_draw(bn1[l], bn2[l], bs1[l], bs2[l], sig1sq, sig2sq, lam2);
        t.nd[l].z = 1;
      }
    }
    for (int i = 0; i < n; ++i) {
      double nv = t.nd[leaf_of[i]].mu;
      E.tot[i] += nv - E.fit[h][i];
      E.fit[h][i] = nv;
    }
  }

  void update_g_hyper() {
    if (G.H == 0) return;
    int K0 = 0, Lg = 0;
    double S0 = 0.0, S1 = 0.0;
    for (int h = 0; h < G.H; ++h) {
      const Tree& t = G.trees[h];
      for (size_t l = 0; l < t.nd.size(); ++l) {
        if (!t.leaf((int)l)) continue;
        ++Lg;
        double th = t.nd[l].mu;
        if (t.nd[l].z == 1) { ++K0; S0 += th * th; } else { S1 += th * th; }
      }
    }
    if (upd_tau0) {
      double nu_post = nu0 + K0;
      double ssq_post = (nu0 * s0sq + S0) / nu_post;
      tau0sq = rsinvchisq_upper(nu_post, ssq_post, tau1sq);
    }
    if (upd_tau1) {
      int K1 = Lg - K0;
      double nu_post = nu1 + K1;
      double ssq_post = (nu1 * s1sq + S1) / nu_post;
      tau1sq = rsinvchisq_lower(nu_post, ssq_post, tau0sq);
    }
    if (!w_fixed) {
      w = R::rbeta(aw + K0, bw + (Lg - K0));
      // guard the exact endpoints so log(w), log(1 - w) stay finite
      if (w < 1e-12) w = 1e-12;
      if (w > 1.0 - 1e-12) w = 1.0 - 1e-12;
    }
  }

  void update_sigma() {
    double ss1 = 0.0, ss2 = 0.0;
    for (int i = 0; i < n; ++i) {
      double mu = F.tot[i] + ((src[i] == 1 && G.H > 0) ? G.tot[i] : 0.0);
      double r = y[i] - mu;
      if (src[i] == 1) ss1 += r * r; else ss2 += r * r;
    }
    if (n2obs > 0)
      sig2sq = (nu_sig * lam_sig2 + ss2) / R::rchisq(nu_sig + n2obs);
    if (n1obs > 0)
      sig1sq = (nu_sig * lam_sig1 + ss1) / R::rchisq(nu_sig + n1obs);
    else
      sig1sq = single_arm ? sig2sq : (nu_sig * lam_sig1) / R::rchisq(nu_sig);
  }

  void impute_censored() {
    if (!aft) return;
    for (int i = 0; i < n; ++i) {
      if (status[i] == 1) continue;
      double mu = F.tot[i] + ((src[i] == 1 && G.H > 0) ? G.tot[i] : 0.0);
      double sg = std::sqrt(src[i] == 1 ? sig1sq : sig2sq);
      y[i] = rtnorm_lower(mu, sg, logC[i]);
    }
  }

  // per-leaf summaries of the g ensemble at the current state:
  // tree, n1, ybar, v1, z, theta, P(z = 1 | ybar)
  NumericMatrix g_leaf_summary() {
    std::vector<double> rows;
    int nl = 0;
    for (int h = 0; h < G.H; ++h) {
      compute_resid(G, h);
      stats(G.trees[h], G);
      const Tree& t = G.trees[h];
      for (size_t l = 0; l < t.nd.size(); ++l) {
        if (!t.leaf((int)l)) continue;
        double n1 = bn1[l], yb = n1 > 0 ? bs1[l] / n1 : NA_REAL;
        rows.push_back(h); rows.push_back(n1); rows.push_back(yb);
        rows.push_back(n1 > 0 ? sig1sq / n1 : NA_REAL);
        rows.push_back(t.nd[l].z); rows.push_back(t.nd[l].mu);
        rows.push_back(g_leaf_pspike(n1, bs1[l], sig1sq, tau0sq, tau1sq, w));
        ++nl;
      }
    }
    NumericMatrix M(nl, 7);
    for (int r = 0; r < nl; ++r) for (int c = 0; c < 7; ++c) M(r, c) = rows[r * 7 + c];
    return M;
  }

  List serialize(const Ensemble& E) const {
    int total = 0;
    for (int h = 0; h < E.H; ++h) total += (int)E.trees[h].nd.size();
    NumericMatrix M(total, 6);
    IntegerVector roots(E.H);
    int off = 0;
    for (int h = 0; h < E.H; ++h) {
      const Tree& t = E.trees[h];
      roots[h] = off;
      for (size_t k = 0; k < t.nd.size(); ++k) {
        const Node& q = t.nd[k];
        M(off + k, 0) = q.var;
        M(off + k, 1) = (q.var >= 0) ? cp[q.var][q.cut] : NA_REAL;
        M(off + k, 2) = (q.left >= 0) ? q.left + off : -1;
        M(off + k, 3) = (q.right >= 0) ? q.right + off : -1;
        M(off + k, 4) = q.mu;
        M(off + k, 5) = q.z;
      }
      off += (int)t.nd.size();
    }
    return List::create(Named("nodes") = M, Named("roots") = roots);
  }
};

void init_ensemble(Ensemble& E, bool isg, int H, double alpha, double beta, int nmin, int n) {
  E.isg = isg; E.H = H; E.alpha = alpha; E.beta = beta; E.nmin = nmin;
  E.trees.assign(H, Tree());
  E.fit.assign(H, std::vector<double>(n, 0.0));
  E.tot.assign(n, 0.0);
  for (int k = 0; k < 3; ++k) { E.acc[k] = 0; E.prop[k] = 0; }
}

CutPoints cutpoints_from_list(List cpl) {
  CutPoints cp(cpl.size());
  for (int v = 0; v < cpl.size(); ++v) {
    NumericVector c = cpl[v];
    cp[v].assign(c.begin(), c.end());
  }
  return cp;
}

// prior tree draw (for c_g Monte Carlo)
void grow_prior(Tree& t, int node, int d, double alpha, double beta,
                const CutPoints& cp, int maxdepth, std::vector<int>& buf) {
  if (d >= maxdepth) return;
  if (unif_rand() >= alpha * std::pow(1.0 + d, -beta)) return;
  avail_vars(t, node, cp, buf);
  if (buf.empty()) return;
  int v = buf[(int)(unif_rand() * buf.size())];
  int lo, hi;
  t.bounds(node, v, (int)cp[v].size(), lo, hi);
  int c = lo + (int)(unif_rand() * (hi - lo + 1));
  t.grow(node, v, c);
  int l = t.nd[node].left, r = t.nd[node].right;
  grow_prior(t, l, d + 1, alpha, beta, cp, maxdepth, buf);
  grow_prior(t, r, d + 1, alpha, beta, cp, maxdepth, buf);
}

// Load an initial forest (warm start) from the serialize() layout, except
// that column 1 holds the cut INDEX into cp[var] rather than the cut value.
// Trees beyond roots.size() stay stumps; fit and tot are set from the leaves.
void load_forest(Ensemble& E, List init, const std::vector<double>& Xv, int n,
                 const CutPoints& cp) {
  NumericMatrix M = init["nodes"];
  IntegerVector roots = init["roots"];
  for (int h = 0; h < E.H && h < roots.size(); ++h) {
    std::vector<int> q; q.push_back(roots[h]);
    std::map<int, int> loc;
    for (size_t k = 0; k < q.size(); ++k) {
      int g = q[k]; loc[g] = (int)k;
      if (M(g, 2) >= 0) { q.push_back((int)M(g, 2)); q.push_back((int)M(g, 3)); }
    }
    Tree t; t.nd.clear();
    for (size_t k = 0; k < q.size(); ++k) {
      int g = q[k]; Node a;
      a.var = (int)M(g, 0); a.cut = (int)M(g, 1);
      a.left = (M(g, 2) >= 0) ? loc[(int)M(g, 2)] : -1;
      a.right = (M(g, 3) >= 0) ? loc[(int)M(g, 3)] : -1;
      a.parent = -1; a.mu = M(g, 4); a.z = (int)M(g, 5);
      if (a.left < 0) { a.var = -1; a.cut = -1; }
      t.nd.push_back(a);
    }
    for (size_t k = 0; k < t.nd.size(); ++k)
      if (t.nd[k].left >= 0) { t.nd[t.nd[k].left].parent = (int)k; t.nd[t.nd[k].right].parent = (int)k; }
    E.trees[h] = t;
    for (int i = 0; i < n; ++i) {
      int l = find_leaf(t, Xv, n, i, cp);
      double v = t.nd[l].mu;
      E.tot[i] += v - E.fit[h][i];
      E.fit[h][i] = v;
    }
  }
}

} // namespace

// ---------------------------------------------------------------------------
// Exported: the sampler
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List lrc_sampler_cpp(NumericMatrix X, NumericVector y, IntegerVector src,
                     IntegerVector status, NumericVector logC, List cutpoints,
                     List hyper, List control) {
  Sampler S;
  S.n = X.nrow(); S.p = X.ncol();
  S.Xv.assign(X.begin(), X.end());
  S.y.assign(y.begin(), y.end());
  S.src.assign(src.begin(), src.end());
  S.status.assign(status.begin(), status.end());
  S.logC.assign(logC.begin(), logC.end());
  S.cp = cutpoints_from_list(cutpoints);
  S.n1obs = 0; S.n2obs = 0;
  for (int i = 0; i < S.n; ++i) { if (S.src[i] == 1) S.n1obs++; else S.n2obs++; }
  S.single_arm = (S.n1obs == 0);

  int Hf = as<int>(hyper["H_f"]), Hg = as<int>(hyper["H_g"]);
  S.lam2 = as<double>(hyper["lambda_f_sq"]);
  S.nu_sig = as<double>(hyper["nu_sigma"]);
  S.lam_sig1 = as<double>(hyper["lambda_sigma1"]);
  S.lam_sig2 = as<double>(hyper["lambda_sigma2"]);
  S.nu0 = as<double>(hyper["nu0"]);
  S.s0sq = as<double>(hyper["s0_sq"]);
  S.nu1 = as<double>(hyper["nu1"]);
  S.s1sq = as<double>(hyper["s1_sq"]);
  S.aw = as<double>(hyper["a_w"]);
  S.bw = as<double>(hyper["b_w"]);
  S.tau1sq = as<double>(hyper["tau1_sq"]);
  S.upd_tau0 = as<bool>(hyper["update_tau0"]);
  S.upd_tau1 = as<bool>(hyper["update_tau1"]);
  S.aft = as<bool>(hyper["aft"]);
  double wfix = as<double>(hyper["w_fixed"]);  // NA to learn
  S.w_fixed = !ISNAN(wfix);
  S.w = S.w_fixed ? wfix : S.aw / (S.aw + S.bw);
  S.tau0sq = std::min(S.s0sq, 0.5 * S.tau1sq);
  S.sig1sq = as<double>(hyper["sigma1_sq_init"]);
  S.sig2sq = as<double>(hyper["sigma2_sq_init"]);
  S.pg = as<double>(hyper["p_grow"]);
  S.pp = as<double>(hyper["p_prune"]);
  S.pc = 1.0 - S.pg - S.pp;
  int nmin_f = as<int>(hyper["n_min"]);
  int nmin_g_req = hyper.containsElementNamed("n_min_g") ? as<int>(hyper["n_min_g"]) : nmin_f;
  int nmin_g = S.single_arm ? 0 : nmin_g_req;

  init_ensemble(S.F, false, Hf, as<double>(hyper["alpha_f"]), as<double>(hyper["beta_f"]), nmin_f, S.n);
  init_ensemble(S.G, true, Hg, as<double>(hyper["alpha_g"]), as<double>(hyper["beta_g"]), nmin_g, S.n);
  S.leaf_of.assign(S.n, 0);
  S.resid.assign(S.n, 0.0);
  // optional warm start of the g forest (cut indices, see load_forest)
  if (Hg > 0 && control.containsElementNamed("g_init") && !Rf_isNull(control["g_init"]))
    load_forest(S.G, control["g_init"], S.Xv, S.n, S.cp);

  int nburn = as<int>(control["n_burn"]), ndraw = as<int>(control["n_draw"]);
  int thin = as<int>(control["thin"]);
  // number of g-ensemble sweeps per iteration (1 = one move per g tree)
  int gsweeps = control.containsElementNamed("g_sweeps") ? as<int>(control["g_sweeps"]) : 1;
  if (gsweeps < 1) gsweeps = 1;
  bool verbose = as<bool>(control["verbose"]);
  bool keep_trees = as<bool>(control["keep_trees"]);
  bool keep_train = as<bool>(control["keep_train"]);
  bool keep_leaf = as<bool>(control["keep_leaf_stats"]);
  List g_leaf(keep_leaf && Hg > 0 ? ndraw : 0);
  int niter = nburn + ndraw * thin;

  List f_draws(keep_trees ? ndraw : 0), g_draws((keep_trees && Hg > 0) ? ndraw : 0);
  NumericVector o_sig1(ndraw), o_sig2(ndraw), o_tau0(ndraw), o_tau1(ndraw), o_w(ndraw);
  IntegerVector o_K0(ndraw), o_Lg(ndraw);
  NumericMatrix f_train(keep_train ? ndraw : 0, S.n), g_train((keep_train && Hg > 0) ? ndraw : 0, S.n);

  // initial imputation for censored observations at the bound plus a little
  if (S.aft) {
    for (int i = 0; i < S.n; ++i) if (S.status[i] == 0) S.y[i] = S.logC[i] + 0.1 * std::sqrt(S.src[i] == 1 ? S.sig1sq : S.sig2sq);
  }

  int d = 0;
  for (int it = 0; it < niter; ++it) {
    for (int h = 0; h < Hf; ++h) S.update_tree(S.F, h);
    for (int s = 0; s < gsweeps; ++s)
      for (int h = 0; h < Hg; ++h) S.update_tree(S.G, h);
    S.update_g_hyper();
    S.update_sigma();
    S.impute_censored();
    if (verbose && ((it + 1) % 100 == 0)) Rprintf("iter %d / %d\n", it + 1, niter);
    if (it >= nburn && ((it - nburn) % thin == 0)) {
      o_sig1[d] = S.sig1sq; o_sig2[d] = S.sig2sq; o_tau0[d] = S.tau0sq; o_tau1[d] = S.tau1sq; o_w[d] = S.w;
      int K0 = 0, Lg = 0;
      for (int h = 0; h < Hg; ++h) {
        const Tree& t = S.G.trees[h];
        for (size_t l = 0; l < t.nd.size(); ++l) if (t.leaf((int)l)) { ++Lg; K0 += t.nd[l].z; }
      }
      o_K0[d] = K0; o_Lg[d] = Lg;
      if (keep_trees) {
        f_draws[d] = S.serialize(S.F);
        if (Hg > 0) g_draws[d] = S.serialize(S.G);
      }
      if (keep_leaf && Hg > 0) g_leaf[d] = S.g_leaf_summary();
      if (keep_train) {
        for (int i = 0; i < S.n; ++i) {
          f_train(d, i) = S.F.tot[i];
          if (Hg > 0) g_train(d, i) = S.G.tot[i];
        }
      }
      ++d;
    }
    if ((it % 50) == 0) Rcpp::checkUserInterrupt();
  }

  NumericMatrix acc(2, 3);
  for (int k = 0; k < 3; ++k) {
    acc(0, k) = S.F.prop[k] > 0 ? (double)S.F.acc[k] / S.F.prop[k] : NA_REAL;
    acc(1, k) = S.G.prop[k] > 0 ? (double)S.G.acc[k] / S.G.prop[k] : NA_REAL;
  }

  return List::create(
    Named("f_draws") = f_draws, Named("g_draws") = g_draws,
    Named("sigma1_sq") = o_sig1, Named("sigma2_sq") = o_sig2,
    Named("tau0_sq") = o_tau0, Named("tau1_sq") = o_tau1, Named("w") = o_w,
    Named("K0") = o_K0, Named("L_g") = o_Lg,
    Named("f_train") = f_train, Named("g_train") = g_train,
    Named("accept") = acc, Named("single_arm") = S.single_arm,
    Named("g_leaf_stats") = g_leaf);
}

// ---------------------------------------------------------------------------
// Exported: prediction from serialized forests
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List predict_forest_cpp(List draws, NumericMatrix Xnew) {
  int nd = draws.size(), m = Xnew.nrow();
  NumericMatrix val(nd, m);
  IntegerMatrix nspike(nd, m);
  for (int d = 0; d < nd; ++d) {
    List dr = draws[d];
    NumericMatrix M = dr["nodes"];
    IntegerVector roots = dr["roots"];
    int H = roots.size();
    for (int i = 0; i < m; ++i) {
      double s = 0.0; int ns = 0;
      for (int h = 0; h < H; ++h) {
        int k = roots[h];
        while (M(k, 2) >= 0) {
          int v = (int)M(k, 0);
          k = (Xnew(i, v) < M(k, 1)) ? (int)M(k, 2) : (int)M(k, 3);
        }
        s += M(k, 4);
        ns += (int)M(k, 5);
      }
      val(d, i) = s; nspike(d, i) = ns;
    }
  }
  return List::create(Named("value") = val, Named("nspike") = nspike);
}

// Estimand-level discrepancy variance per draw:
//   V_g = sum_h sum_{leaves l of h} (n_l / N)^2 tau^2_{z_hl},
// with n_l the number of the N profiles in leaf l (the per-leaf form of
// c_g sum_h tau^2_{z_h} in spec Section C).
// [[Rcpp::export]]
NumericVector forest_gvar_cpp(List draws, NumericMatrix Xnew, NumericVector tau0sq,
                              NumericVector tau1sq) {
  int nd = draws.size(), m = Xnew.nrow();
  NumericVector out(nd);
  for (int d = 0; d < nd; ++d) {
    List dr = draws[d];
    NumericMatrix M = dr["nodes"];
    IntegerVector roots = dr["roots"];
    std::vector<double> cnt(M.nrow(), 0.0);
    for (int i = 0; i < m; ++i) {
      for (int h = 0; h < roots.size(); ++h) {
        int k = roots[h];
        while (M(k, 2) >= 0) {
          int v = (int)M(k, 0);
          k = (Xnew(i, v) < M(k, 1)) ? (int)M(k, 2) : (int)M(k, 3);
        }
        cnt[k] += 1.0;
      }
    }
    double V = 0.0;
    for (int k = 0; k < M.nrow(); ++k) {
      if (M(k, 2) >= 0 || cnt[k] == 0.0) continue;
      double sh = cnt[k] / m;
      V += sh * sh * (M(k, 5) > 0.5 ? tau0sq[d] : tau1sq[d]);
    }
    out[d] = V;
  }
  return out;
}

// ---------------------------------------------------------------------------
// Exported: c_g Monte Carlo under the g-tree prior
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double cg_montecarlo_cpp(NumericMatrix X, List cutpoints, double alpha, double beta,
                         int nsim, int maxdepth) {
  CutPoints cp = cutpoints_from_list(cutpoints);
  int n = X.nrow();
  std::vector<double> Xv(X.begin(), X.end());
  std::vector<int> buf;
  double acc = 0.0;
  for (int s = 0; s < nsim; ++s) {
    Tree t;
    grow_prior(t, 0, 0, alpha, beta, cp, maxdepth, buf);
    std::vector<double> cnt(t.nd.size(), 0.0);
    for (int i = 0; i < n; ++i) cnt[find_leaf(t, Xv, n, i, cp)] += 1.0;
    double c = 0.0;
    for (size_t l = 0; l < cnt.size(); ++l) c += (cnt[l] / n) * (cnt[l] / n);
    acc += c;
  }
  return acc / nsim;
}

// ---------------------------------------------------------------------------
// Exported: unit-test helpers for the leaf-level pieces
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double f_leaf_logmarg_cpp(double n1, double n2, double s1, double s2,
                          double sig1sq, double sig2sq, double lam2) {
  return f_leaf_lm(n1, n2, s1, s2, sig1sq, sig2sq, lam2);
}

// [[Rcpp::export]]
double g_leaf_logmarg_cpp(double n1, double s1, double sig1sq, double tau0sq,
                          double tau1sq, double w) {
  return g_leaf_lm(n1, s1, sig1sq, tau0sq, tau1sq, w);
}

// [[Rcpp::export]]
double g_leaf_pspike_cpp(double n1, double s1, double sig1sq, double tau0sq,
                         double tau1sq, double w) {
  return g_leaf_pspike(n1, s1, sig1sq, tau0sq, tau1sq, w);
}

// [[Rcpp::export]]
List g_leaf_draw_cpp(double n1, double s1, double sig1sq, double tau0sq,
                     double tau1sq, double w, int ndraw) {
  IntegerVector z(ndraw); NumericVector th(ndraw);
  for (int d = 0; d < ndraw; ++d) {
    int zz; double tt;
    g_leaf_draw(n1, s1, sig1sq, tau0sq, tau1sq, w, zz, tt);
    z[d] = zz; th[d] = tt;
  }
  return List::create(Named("z") = z, Named("theta") = th);
}

// [[Rcpp::export]]
NumericVector update_tau0_cpp(NumericVector theta, IntegerVector z, double nu0,
                              double s0sq, double tau1sq, int ndraw) {
  int K0 = 0; double S0 = 0.0;
  for (int i = 0; i < theta.size(); ++i) if (z[i] == 1) { ++K0; S0 += theta[i] * theta[i]; }
  NumericVector out(ndraw);
  double nu_post = nu0 + K0, ssq_post = (nu0 * s0sq + S0) / nu_post;
  for (int d = 0; d < ndraw; ++d) out[d] = rsinvchisq_upper(nu_post, ssq_post, tau1sq);
  return out;
}

// [[Rcpp::export]]
NumericVector update_tau1_cpp(NumericVector theta, IntegerVector z, double nu1,
                              double s1sq, double tau0sq, int ndraw) {
  int K1 = 0; double S1 = 0.0;
  for (int i = 0; i < theta.size(); ++i) if (z[i] == 0) { ++K1; S1 += theta[i] * theta[i]; }
  NumericVector out(ndraw);
  double nu_post = nu1 + K1, ssq_post = (nu1 * s1sq + S1) / nu_post;
  for (int d = 0; d < ndraw; ++d) out[d] = rsinvchisq_lower(nu_post, ssq_post, tau0sq);
  return out;
}

// [[Rcpp::export]]
NumericVector update_w_cpp(IntegerVector z, double aw, double bw, int ndraw) {
  int K0 = 0, L = z.size();
  for (int i = 0; i < L; ++i) K0 += z[i];
  NumericVector out(ndraw);
  for (int d = 0; d < ndraw; ++d) out[d] = R::rbeta(aw + K0, bw + (L - K0));
  return out;
}

// [[Rcpp::export]]
NumericVector rtnorm_lower_cpp(double mu, double sigma, double a, int n) {
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = rtnorm_lower(mu, sigma, a);
  return out;
}
