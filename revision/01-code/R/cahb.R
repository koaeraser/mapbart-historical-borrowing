# cahb.R: the covariate-adaptive historical control borrowing method (CAHB)
# of Jin, Kim, Scheffler and Jiang (2023), Statistics in Medicine 42(29):
# 5338-5352, doi 10.1002/sim.9913, as a comparator with the harness interface
# of comparators.R: fit_cahb(dat, outcome, ...) returns make_fit()'s list.
#
# Source of the formulation. The paper (Sections 2 and 4, supplement A.2)
# and the authors' code, github.com/JINhuaqing/HistTrial (utils.R,
# simplex.R, realData/linearSimuRealData*.R), which we read in full. No
# package exists (CRAN and GitHub searched); the repository is a set of
# simulation scripts, so the estimator is re-implemented here. Only the
# analysis part of CAHB is relevant: the design part (kernel biased-coin
# allocation driven by the covariate-dependent effective control sample size)
# has no role in a completed RCT with fixed allocation.
#
# The model (paper Section 2, normal likelihood). Current-trial data (Y_i,
# X_i, Z_i), i = 1..n; theta_0(X) the historical control mean function;
#   mu_0(X) ~ N(theta_0(X), 1 / tau(X)),  tau(X)^{1/2} ~ half-normal(gamma),
#   phi_0^2, phi_1^2 ~ Inv-Gamma(0.01, 0.01),  flat prior on mu_1(X).
# Kernel-local posterior modes (Section 4, eq. 8-13), iterated to
# convergence (Algorithm 1), with kernel weights w_ij = K_H(X_i - X_j):
#   mu_0(X_j) = sum_i (1-Z_i) w_ij {Y_i / phi_0^2 + tau_i theta_0(X_i)} /
#               sum_i (1-Z_i) w_ij {1 / phi_0^2 + tau_i}
#   tau~(X_j) = sum_i (1-Z_i) w_ij / [sum_i (1-Z_i) w_ij {mu_0(X_i) -
#               theta_0(X_i)}^2 + 1 / gamma^2]
#   tau       = Proj{ tau~ 1(tau~ > lambda_1); lambda_2 }, the Euclidean
#               projection of the thresholded vector onto the L1 ball of
#               radius lambda_2 (Duchi et al. 2008), as in the authors' code
#   phi_0^2   = {sum_i (1-Z_i) (Y_i - mu_0(X_i))^2 / 2 + b} / (a + n_0/2 + 1)
#   mu_1(X_j) = sum_i Z_i w_ij Y_i / sum_i Z_i w_ij,  phi_1^2 likewise.
# Posterior of mu_0(X) given the plug-in estimates (supplement A.2):
#   N(mu_0(X), 1 / A(X)),  A(X) = sum_i (1-Z_i) w_ij {1/phi_0^2 + tau_i},
# and the information ratio R(X) = A(X) phi_0^2 / sum_i (1-Z_i) w_ij, whose
# excess R(X) - 1 is the historical share of the local control information.
# The covariate-dependent effective control sample size (CECSS) is R(X)
# times the local control count.
#
# Hyperparameters as published: gamma = sqrt(3); lambda_2 = 300 log(n);
# lambda_1 = the 10% quantile of tau~ under self-consistency (Algorithm 2:
# theta_0 replaced by the no-borrowing kernel fit of the current controls,
# lambda_1 = 0); Gaussian product kernel with the diagonal bandwidth matrix
# from Scott's rule of thumb for continuous covariates; a = b = 0.01;
# convergence when the largest mean squared change of (mu_0, tau, phi_0^2)
# is below 1e-5, at most 100 iterations.
#
# Points the paper leaves unspecified, or that do not transfer to this
# study, and the choice made for the primary method "CAHB". Two variants are
# registered alongside: "CAHB-p10" (all ten covariates in the kernel,
# otherwise identical) and "CAHB-lit" (every published constant taken
# literally: linear theta_0, all covariates, Scott's rule, gamma^2 = 3,
# lambda_2 = 300 log n, independent per-profile draws). The choices below
# were settled on three replicates of Sc1, Sc5_d2 and Sc4_X5_d2 (seeds 2027,
# 2050, 2101, none of which are study replicates) by looking at whether the
# estimator operates at all (local sample size, phi_0 estimate, the tau
# switch under a global shift); they were not tuned to the study's metrics.
#  1. theta_0(X). The paper's historical model is a subgroup-wise linear
#     regression fitted to the historical patients (4 covariates, two
#     binary). Here theta_0 is the same kernel regression the method uses for
#     mu_0 (Nadaraya-Watson on the kernel covariates with the same
#     bandwidth), fitted to the RWD controls and evaluated at the RCT
#     profiles. A sharper historical estimate (BART on the RWD) was tried
#     and rejected: paired with the kernel-smoothed mu_0 the difference
#     mu_0 - theta_0 reflects the smoother's bias rather than the shift, so
#     tau is small everywhere and nothing is borrowed even under Sc1. The
#     historical estimate enters as a known function; the paper does not
#     account for its variance (Section 6 states this as a limitation) and
#     neither do we.
#  2. Kernel covariates. The paper matches on four literature-chosen
#     prognostic covariates. The primary method uses the four covariates
#     with the largest variable-inclusion counts in a BART fit to the RWD
#     controls (BART_OPTS of comparators.R; on the diagnostic replicates
#     these were X1, X5, X6 and one of X2 to X4). CAHB-p10 uses all ten.
#  3. Bandwidth. Scott's rule h_k = s_k n^{-1/(p+4)} is the paper's rule for
#     continuous covariates. With p = 10 and the n = 100 current controls it
#     gives h = 0.72 s_k, at which the product kernel localizes to single
#     subjects (median local effective control sample size sum_i w_ij about
#     1.0, of which 1 is the subject itself): the estimator interpolates,
#     phi_0 collapses to about 0.1, profiles without a neighbour get a
#     near-infinite posterior variance, and the tau threshold cannot
#     separate a shift from noise (separation needs sum_i w_ij d^2 much
#     larger than 1/gamma^2). LOO cross-validation of the no-borrowing
#     smoother does not help: in 4 or 10 dimensions it selects 1 to 1.5
#     times Scott's rule (local sample 1 to 5) because Nadaraya-Watson is
#     poor at every bandwidth for this surface, and at that bandwidth CAHB
#     borrowed from the shifted RWD in Sc5_d2 on two of three replicates.
#     The primary method therefore multiplies Scott's rule by the constant
#     at which the median local effective control sample size over the RCT
#     control profiles equals 20 (bisection; about 2.5 in four dimensions,
#     about 3.4 in ten). The kernel is normalized to K(0) = 1 (paper Section
#     3), so sum_i w_ij is a local sample size and phi_0^2 / sum_i w_ij is on
#     the right scale; the authors' code uses an unnormalized Gaussian
#     density, which rescales the posterior variance by (2 pi)^{-p/2}
#     |H|^{-1/2} (about 0.7 in their real-data setting) and changes tau~
#     through the 1/gamma^2 term. The same H serves mu_0, tau and mu_1, as in
#     the paper.
#  4. gamma and lambda_2 are on the outcome scale. The published values are
#     tied to the paper's residual variance phi_0^2 = 0.11 (24-month hip
#     BMD); the authors' own linear-model script (linearSimu.R, phi_0 = 3)
#     uses gamma^2 = 0.01 and lambda_2 = 30 log n instead. The dimensionless
#     quantities are gamma^2 phi_0^2 (the prior-scale cap on the borrowed
#     precision relative to the data precision) and lambda_2 phi_0^2 / n (the
#     projection cap on the mean borrowed-to-data precision ratio). The
#     primary method keeps the published dimensionless values: gamma^2 =
#     3 x 0.11 / phi_hat^2 and lambda_2 = 300 log(n) x 0.11 / phi_hat^2,
#     where phi_hat^2 is the no-borrowing residual variance of the RCT
#     controls under the method's own smoother (the reference model of
#     Algorithm 2; about 4 to 5 at the local-sample-20 bandwidth, larger
#     than the true 2.25 because the smoother leaves structure in the
#     residual), so gamma^2 is about 0.07 and lambda_2 about 40.
#  5. Posterior draws of the averaged estimand. The authors' code draws
#     mu_0(X_j) and mu_1(X_j) independently across profiles j from their
#     marginal posteriors and averages; the resulting posterior of the
#     profile-averaged difference is not calibrated (the paper calibrates
#     the test cut-off c by simulation instead of relying on it). The
#     primary method keeps every marginal posterior of the paper (mean
#     mu_0(X_j), variance 1/A(X_j)) and induces the dependence across
#     profiles through the kernel smoother's own weights: with a_ij =
#     (1-Z_i) w_ij / phi_0^2 / A(X_j) and standard normal epsilon_i, one draw
#     is mu_0(X_j) + A(X_j)^{-1/2} sum_i a_ij epsilon_i / (sum_i a_ij^2)^{1/2};
#     likewise for mu_1. Point estimates, the surface and the map are
#     unaffected; only the posterior SD, coverage and power of the averaged
#     estimands depend on this.
#  6. Treated arm and standardization: mu_1 from the RCT treated subjects
#     only (flat prior), the estimands standardized over all n RCT profiles
#     as for every other method (the paper's point estimate delta_hat =
#     N^{-1} sum_i {mu_1(X_i) - mu_0(X_i)} is the same quantity).
#  7. tau is evaluated at all n current-trial profiles (treated and control),
#     and the L1 projection acts on that length-n vector with radius
#     lambda_2, as in the authors' mOptTau (radius lam * log(m), m = n).
#  8. Harness outputs: map = 1 - 1/R(X) at the RCT control profiles (the
#     share of the local posterior precision contributed by the historical
#     prior, in [0, 1)); ess = list(realized = sum over RCT controls of
#     R(X_i) - 1, the CECSS excess over the actual control count; prior is
#     NULL because CAHB has no data-free prior ESS); extra holds tau, R,
#     lambda_1, lambda_2, gamma^2, phi_0, the bandwidth, the kernel
#     covariates and their BART counts, the local sample sizes and the
#     iteration count.
#
# Gaussian outcome only (the paper treats a continuous outcome; the
# survival study is not attempted).

# ---------------------------------------------------------------------------
# Kernel helpers
# ---------------------------------------------------------------------------

# Gaussian product kernel with K(0) = 1: n x m matrix of K_h(X_i - Q_j).
cahb_kernel <- function(X, Q, h) {
  Xs <- sweep(X, 2, h, "/"); Qs <- sweep(Q, 2, h, "/")
  d2 <- outer(rowSums(Xs^2), rowSums(Qs^2), "+") - 2 * tcrossprod(Xs, Qs)
  exp(-0.5 * pmax(d2, 0))
}

# Scott's rule of thumb: h_k = s_k n^{-1/(p+4)}.
cahb_scott_bw <- function(X) apply(X, 2, sd) * nrow(X)^(-1 / (ncol(X) + 4))

# Nadaraya-Watson fit of y on X evaluated at Q (with a nearest-neighbour
# fallback where the kernel mass underflows).
cahb_nw <- function(y, X, Q, h) {
  W <- cahb_kernel(X, Q, h)
  den <- colSums(W)
  est <- as.vector(crossprod(W, y)) / den
  bad <- !is.finite(est) | den < 1e-300
  if (any(bad)) {
    Xs <- sweep(X, 2, h, "/"); Qs <- sweep(Q[bad, , drop = FALSE], 2, h, "/")
    d2 <- outer(rowSums(Qs^2), rowSums(Xs^2), "+") - 2 * tcrossprod(Qs, Xs)
    est[bad] <- y[apply(d2, 1, which.min)]
  }
  est
}

# Leave-one-out cross-validated multiplier of the rule-of-thumb bandwidth
# for the Nadaraya-Watson smoother of y on X.
cahb_cv_mult <- function(y, X, h0, grid = c(1, 1.5, 2, 3, 4, 6, 8, 12)) {
  sse <- vapply(grid, function(cc) {
    W <- cahb_kernel(X, X, cc * h0)
    diag(W) <- 0
    den <- colSums(W)
    pred <- as.vector(crossprod(W, y)) / den
    pred[!is.finite(pred) | den < 1e-300] <- mean(y)
    sum((y - pred)^2)
  }, 0)
  grid[which.min(sse)]
}

# Euclidean projection onto the L1 ball of radius s (Duchi et al. 2008), as
# in the authors' simplex.R. Entries are nonnegative here.
cahb_proj_l1 <- function(v, s) {
  if (sum(abs(v)) <= s) return(v)
  u <- sort(abs(v), decreasing = TRUE)
  css <- cumsum(u)
  rho <- which(u * seq_along(u) > (css - s))
  if (length(rho) == 0) return(rep(0, length(v)))
  rho <- rho[length(rho)]
  theta <- (css[rho] - s) / rho
  w <- pmax(abs(v) - theta, 0)
  w * sign(v)
}

# ---------------------------------------------------------------------------
# Algorithm 1: kernel-local posterior modes of (mu_0, tau, phi_0)
# ---------------------------------------------------------------------------
# W: n x n kernel matrix (rows: data points i, columns: query points j);
# y, Z: length n; theta0: theta_0 at the n current-trial profiles.
# Returns mu0 (length n), tau (length n, thresholded and projected),
# tau_raw (before threshold/projection), phi0, A, R, and the iteration count.
cahb_alg1 <- function(y, Z, W, theta0, gamma2, lam2, lam1, a0 = 0.01, b0 = 0.01,
                      maxit = 100, tol = 1e-5, phi0_init = 1) {
  n <- length(y); ctrl <- Z == 0; n0 <- sum(ctrl)
  Wc <- W[ctrl, , drop = FALSE]
  yc <- y[ctrl]; th_c <- theta0[ctrl]
  sw <- colSums(Wc)
  tau <- rep(0, n); phi0 <- phi0_init
  mu0 <- NULL; it <- 0
  repeat {
    it <- it + 1
    tc <- tau[ctrl]
    A <- colSums(Wc * (1 / phi0^2 + tc))
    Bv <- colSums(Wc * (yc / phi0^2 + tc * th_c))
    mu0_new <- Bv / A
    phi0_new <- sqrt((sum((yc - mu0_new[ctrl])^2) / 2 + b0) / (1 + a0 + n0 / 2))
    d2 <- (mu0_new[ctrl] - th_c)^2
    tau_raw <- sw / (colSums(Wc * d2) + 1 / gamma2)
    tau_new <- tau_raw
    tau_new[tau_new <= lam1] <- 0
    tau_new <- cahb_proj_l1(tau_new, lam2)
    if (!is.null(mu0)) {
      err <- max(mean((mu0_new - mu0)^2), mean((tau_new - tau)^2), (phi0_new^2 - phi0^2)^2)
      mu0 <- mu0_new; tau <- tau_new; phi0 <- phi0_new
      if (err < tol || it >= maxit) break
    } else {
      mu0 <- mu0_new; tau <- tau_new; phi0 <- phi0_new
    }
  }
  tc <- tau[ctrl]
  A <- colSums(Wc * (1 / phi0^2 + tc))
  R <- A * phi0^2 / sw
  list(mu0 = mu0, tau = tau, tau_raw = tau_raw, phi0 = phi0, A = A, R = R, iter = it)
}

# Draws of a kernel-smoothed surface: joint (shared latent residuals, the
# paper's marginals) or independent per profile (the authors' code).
cahb_draws <- function(m, v, Wsub, weights, B, joint) {
  n <- length(m)
  if (!joint) {
    return(matrix(m, B, n, byrow = TRUE) + matrix(rnorm(B * n), B, n) * rep(sqrt(v), each = B))
  }
  # a_ij proportional to Wsub_ij * weights_i (rows i: source subjects).
  a <- Wsub * weights
  a <- sweep(a, 2, sqrt(colSums(a^2)), "/")
  E <- matrix(rnorm(B * nrow(a)), B, nrow(a))
  matrix(m, B, n, byrow = TRUE) + (E %*% a) * rep(sqrt(v), each = B)
}

# ---------------------------------------------------------------------------
# Bandwidth rule by local effective sample size, covariate ranking, theta_0
# ---------------------------------------------------------------------------
# Multiplier c of the rule-of-thumb bandwidth such that the median over the
# profiles Q of the local effective sample size sum_i K_{c h0}(X_i - Q_j)
# equals target (monotone in c; bisection on log c).
cahb_ess_mult <- function(X, Q, h0, target = 20, lo = 0.25, hi = 50) {
  f <- function(lc) median(colSums(cahb_kernel(X, Q, exp(lc) * h0))) - target
  if (f(log(lo)) >= 0) return(lo)
  if (f(log(hi)) <= 0) return(hi)
  exp(uniroot(f, c(log(lo), log(hi)), tol = 1e-3)$root)
}

# Covariate ranking by BART variable-inclusion counts in a fit to the RWD
# controls (BART_OPTS of comparators.R).
cahb_rank_bart <- function(dat) {
  Xw <- xmat(dat$rwd_ctrl); yw <- dat$rwd_ctrl$y
  fit <- capture_bart(wbart(x.train = Xw, y.train = yw,
                            ntree = BART_OPTS$ntree, k = BART_OPTS$k, power = BART_OPTS$power,
                            base = BART_OPTS$base, sigdf = BART_OPTS$sigdf,
                            sigquant = BART_OPTS$sigquant, ndpost = BART_OPTS$ndpost,
                            nskip = BART_OPTS$nskip, printevery = 100000L))
  colMeans(fit$varcount)
}

# Historical mean function theta_0 at the RCT profiles X. "nw": kernel
# regression of the RWD control outcomes on the kernel covariates vars with
# the same bandwidth h as the current-trial smoother; "lm": linear
# regression on all covariates (the paper's own historical model).
cahb_theta0 <- function(dat, X, method, vars, h) {
  Xw <- xmat(dat$rwd_ctrl); yw <- dat$rwd_ctrl$y
  if (method == "nw") return(cahb_nw(yw, Xw[, vars, drop = FALSE], X[, vars, drop = FALSE], h))
  if (method == "lm") {
    beta <- coef(lm.fit(cbind(1, Xw), yw))
    return(as.vector(cbind(1, X) %*% beta))
  }
  stop("unknown theta0_method: ", method)
}

# ---------------------------------------------------------------------------
# The harness method
# ---------------------------------------------------------------------------
# kernel_vars: "top" uses the n_top covariates with the largest BART
# variable-inclusion counts in the RWD fit (the paper matches on four
# literature-chosen prognostic covariates), "all" uses every covariate.
# bw_rule: "ess" multiplies Scott's rule so that the median local effective
# control sample size at the RCT control profiles equals ess_target; "cv"
# uses the LOO-CV multiplier of the no-borrowing smoother; "scott" uses the
# rule of thumb as is.
fit_cahb <- function(dat, outcome = "gaussian", theta0_method = c("nw", "lm"),
                     kernel_vars = c("top", "all"), n_top = 4,
                     bw_rule = c("ess", "cv", "scott"), ess_target = 20,
                     gamma2_pub = 3, lam2_coef = 300, phi2_ref = 0.11, lam1_q = 0.10,
                     bw_grid = c(1, 1.5, 2, 3, 4, 6, 8, 12),
                     joint_draws = TRUE, n_draw = N_DRAW, ...) {
  if (outcome != "gaussian") stop("CAHB is defined for the Gaussian study only")
  theta0_method <- match.arg(theta0_method); kernel_vars <- match.arg(kernel_vars)
  bw_rule <- match.arg(bw_rule)
  X <- dat$x_rct; Z <- dat$z_rct
  y <- c(dat$rct_trt$y, dat$rct_ctrl$y)    # same order as x_rct (treated first)
  stopifnot(length(y) == nrow(X), all(Z == rep(c(1L, 0L), c(nrow(dat$rct_trt), nrow(dat$rct_ctrl)))))
  n <- nrow(X); ctrl <- Z == 0; trt <- !ctrl

  # 1. Kernel covariates (deviation 2) and bandwidth (deviation 3).
  varcount <- NULL
  if (kernel_vars == "top") {
    varcount <- cahb_rank_bart(dat)
    vars <- sort(order(varcount, decreasing = TRUE)[seq_len(min(n_top, ncol(X)))])
  } else {
    vars <- seq_len(ncol(X))
  }
  Xk <- X[, vars, drop = FALSE]
  h0 <- cahb_scott_bw(Xk)
  mult_c <- switch(bw_rule,
                   ess = cahb_ess_mult(Xk[ctrl, , drop = FALSE], Xk[ctrl, , drop = FALSE], h0, ess_target),
                   cv = cahb_cv_mult(y[ctrl], Xk[ctrl, , drop = FALSE], h0, bw_grid),
                   scott = 1)
  h <- mult_c * h0
  W <- cahb_kernel(Xk, Xk, h)

  # 2. Historical mean function theta_0 at the RCT profiles (deviation 1).
  theta0 <- cahb_theta0(dat, X, theta0_method, vars, h)

  # 3. Reference (no-borrowing) fit: phi_hat^2 and the self-consistent
  #    theta_0 for Algorithm 2.
  ref <- cahb_alg1(y, Z, W, theta0 = rep(0, n), gamma2 = 1, lam2 = 0, lam1 = Inf)
  phi2_hat <- ref$phi0^2
  scale <- if (is.null(phi2_ref)) 1 else phi2_ref / phi2_hat
  gamma2 <- gamma2_pub * scale
  lam2 <- lam2_coef * log(n) * scale

  # 4. Algorithm 2: lambda_1 = the lam1_q quantile of tau~ when theta_0 is the
  #    no-borrowing fit of the current controls.
  self <- cahb_alg1(y, Z, W, theta0 = ref$mu0, gamma2 = gamma2, lam2 = lam2, lam1 = 0,
                    phi0_init = ref$phi0)
  lam1 <- unname(quantile(self$tau_raw, lam1_q))

  # 5. Algorithm 1 with the historical theta_0.
  fit0 <- cahb_alg1(y, Z, W, theta0 = theta0, gamma2 = gamma2, lam2 = lam2, lam1 = lam1,
                    phi0_init = ref$phi0)

  # 6. Treated arm: flat-prior kernel fit and phi_1 (eq. 12-13).
  Wt <- W[trt, , drop = FALSE]
  swt <- colSums(Wt)
  mu1 <- as.vector(crossprod(Wt, y[trt])) / swt
  phi1 <- sqrt((sum((y[trt] - mu1[trt])^2) / 2 + 0.01) / (1 + 0.01 + sum(trt) / 2))

  # 7. Posterior draws (deviation 5).
  v0 <- 1 / fit0$A
  v1 <- phi1^2 / swt
  mu0_draws <- cahb_draws(fit0$mu0, v0, W[ctrl, , drop = FALSE],
                          rep(1 / fit0$phi0^2, sum(ctrl)), n_draw, joint_draws)
  mu1_draws <- cahb_draws(mu1, v1, Wt, rep(1, sum(trt)), n_draw, joint_draws)

  Rc <- fit0$R[ctrl]
  make_fit(mu1_draws, mu0_draws, "gaussian",
           map = pmax(0, 1 - 1 / Rc),
           ess = list(prior = NULL, realized = sum(Rc - 1), R_mean = mean(Rc)),
           extra = list(tau = fit0$tau, tau_raw = fit0$tau_raw, R = fit0$R,
                        lam1 = lam1, lam2 = lam2, gamma2 = gamma2, phi0 = fit0$phi0,
                        phi1 = phi1, phi2_hat = phi2_hat, mult_c = mult_c, h = h,
                        vars = vars, varcount = varcount, iter = fit0$iter,
                        theta0 = theta0, mu0_hat = fit0$mu0, mu1_hat = mu1,
                        sum_w_ctrl = colSums(W[ctrl, , drop = FALSE])))
}

# All ten covariates in the kernel, otherwise as the primary method.
fit_cahb_p10 <- function(dat, outcome = "gaussian", ...) {
  fit_cahb(dat, outcome, kernel_vars = "all")
}

# Literal transfer of the published constants (see the header): linear
# theta_0, all covariates, Scott's rule, gamma^2 = 3, lambda_2 = 300 log n,
# independent per-profile draws.
fit_cahb_lit <- function(dat, outcome = "gaussian", ...) {
  fit_cahb(dat, outcome, theta0_method = "lm", kernel_vars = "all", bw_rule = "scott",
           gamma2_pub = 3, lam2_coef = 300, phi2_ref = NULL, joint_draws = FALSE)
}

if (exists("register_method")) {
  register_method("CAHB", fit_cahb)
  register_method("CAHB-p10", fit_cahb_p10)
  register_method("CAHB-lit", fit_cahb_lit)
}
