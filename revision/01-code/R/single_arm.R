# single_arm.R: the paper's single-arm design (Section "Single-Arm Trial
# Borrowing"): the RCT contributes a treatment arm only, the RWD (n = 300) is
# the sole control source, and the control surface at the trial profiles is
# reconstructed from the RWD. Sourced after harness.R / source_all() and
# fit_lrcbart.R; nothing in the existing files is changed.
#
# Data. gen_data() is reused with n_rct = n_trt and a randomization ratio so
# large that ratio / (ratio + 1) equals 1 in floating point, so every RCT
# subject is treated (z = 1) and rct_ctrl has zero rows; the RNG stream and
# the seed pairing of the two-arm harness are unchanged. single_arm_gen_args()
# returns the gen_args list for run_sim().
#
# Methods (registered by register_single_arm_methods()):
#   LM-CP / AFT-CP and BART-CP: the existing fitters, which already handle an
#     empty RCT control set (the pooled control data are the RWD alone).
#   HierLM / HierAFT: fit_hier_single(), the existing hierarchical model with
#     an empty group 1: theta_1 | mu, tau^2 ~ N(mu, diag(tau^2)) has no data
#     and is drawn from the hierarchy (the paper's synthetic-control reading),
#     and the RCT residual scale, which has no data, is set to the RWD
#     residual scale (sigma_1 = sigma_2 as in the paper's application). With
#     RCT controls present it delegates to fit_hier().
#   LRC-BART in single-arm mode (lrcbart's n_1 = 0 mode: f is BART on the RWD,
#     g a prior draw, spec Section A1), w fixed at 1 (the reported model) or
#     0.9 (sensitivity), N_target at f90 or the absolute 100 (capped at the
#     ESS ceiling when infeasible), and LRC-BART-full: w = 1 with s_0^2 at the
#     smallest grid value 1e-6, the paper's "near-complete pooling".
#   The compatibility map is empty in single-arm mode (no RCT control
#   profiles), so the map columns are NA; ESS_0, the realized ESS, the ceiling
#   and the cap flag are reported as in the two-arm study.

single_arm_gen_args <- function(n_trt, n_rwd = 300) {
  list(n_rct = n_trt, n_rwd = n_rwd, ratio = .Machine$double.xmax)
}

# Hierarchical model with an empty RCT control group (see header).
gibbs_hier_single <- function(y2, X2, status2 = NULL, nu = 3, q = 0.90, hyper_sd = 10,
                              tau_a = 1.5, tau_b = 0.75, n_burn = N_BURN, n_draw = N_DRAW) {
  K <- ncol(X2) + 1
  Xd <- cbind(1, X2); n2 <- length(y2)
  cens <- if (is.null(status2)) rep(FALSE, n2) else status2 == 0
  y_aug <- y2
  XtX <- crossprod(Xd)
  sig_hat <- sigma_prelim(y2, X2, status2)
  lambda <- lambda_bart(sig_hat, nu, q)
  theta2 <- qr.solve(Xd, y2)
  sigma2 <- sig_hat^2
  mu <- theta2
  tau2 <- rep(0.5, K)
  theta1 <- mu
  out_theta <- matrix(NA_real_, n_draw, K)
  out_sigma <- numeric(n_draw)
  for (it in seq_len(n_burn + n_draw)) {
    if (any(cens)) {
      mu_c <- as.vector(Xd[cens, , drop = FALSE] %*% theta2)
      y_aug[cens] <- rtnorm_lower(mu_c, sqrt(sigma2), y2[cens])
    }
    V_inv <- XtX / sigma2 + diag(1 / tau2)
    R <- chol(V_inv)
    rhs <- crossprod(Xd, y_aug) / sigma2 + mu / tau2
    m <- backsolve(R, backsolve(R, rhs, transpose = TRUE))
    theta2 <- as.vector(m + backsolve(R, rnorm(K)))
    ssr <- sum((y_aug - Xd %*% theta2)^2)
    sigma2 <- rinvgamma1((nu + n2) / 2, (nu * lambda + ssr) / 2)
    # group 1 has no data: its coefficients are a draw from the hierarchy
    theta1 <- rnorm(K, mu, sqrt(tau2))
    prec <- 2 / tau2 + 1 / hyper_sd^2
    mean_mu <- ((theta1 + theta2) / tau2) / prec
    mu <- rnorm(K, mean_mu, sqrt(1 / prec))
    ss <- (theta1 - mu)^2 + (theta2 - mu)^2
    tau2 <- 1 / rgamma(K, shape = tau_a + 1, rate = tau_b + ss / 2)
    if (it > n_burn) {
      out_theta[it - n_burn, ] <- theta1
      out_sigma[it - n_burn] <- sqrt(sigma2)
    }
  }
  list(theta = out_theta, sigma = out_sigma)
}

fit_hier_single <- function(dat, outcome, ...) {
  if (nrow(dat$rct_ctrl) > 0) return(fit_hier(dat, outcome, ...))
  trt <- fit_trt_lm(dat, outcome)
  d2 <- prep_xy(dat$rwd_ctrl, outcome)
  f <- gibbs_hier_single(d2$y, d2$X, d2$status)
  mu0 <- predict_lm(f, dat$x_rct)
  make_fit(trt$mu, mu0, outcome, trt$sigma, f$sigma, dat$constants$tau)
}

# LRC-BART single-arm variants (fit_lrcbart handles n_1 = 0; the sampler
# checks it through single_arm = TRUE).
fit_lrcbart_sa <- function(dat, outcome, N_target, w, calib_grid = NULL, ...) {
  if (nrow(dat$rct_ctrl) > 0) stop("single-arm LRC-BART called with RCT controls present")
  fit_lrcbart(dat, outcome, N_target = N_target, w_fixed = w, calib_grid = calib_grid,
              single_arm = TRUE, ...)
}
fit_lrcbart_sa_f90_w1 <- function(dat, outcome, ...) fit_lrcbart_sa(dat, outcome, c("frac", 0.9), 1, ...)
fit_lrcbart_sa_100_w1 <- function(dat, outcome, ...) fit_lrcbart_sa(dat, outcome, 100, 1, ...)
fit_lrcbart_sa_f90_w09 <- function(dat, outcome, ...) fit_lrcbart_sa(dat, outcome, c("frac", 0.9), 0.9, ...)
fit_lrcbart_sa_100_w09 <- function(dat, outcome, ...) fit_lrcbart_sa(dat, outcome, 100, 0.9, ...)
# near-complete pooling: s_0^2 at the smallest grid value (the point-mass limit)
fit_lrcbart_sa_full <- function(dat, outcome, ...) fit_lrcbart_sa(dat, outcome, 100, 1, calib_grid = 1e-6, ...)

register_single_arm_methods <- function() {
  register_method("HierLM", fit_hier_single)
  register_method("HierAFT", fit_hier_single)
  register_method("LRC-BART-f90-w1", fit_lrcbart_sa_f90_w1)
  register_method("LRC-BART-100-w1", fit_lrcbart_sa_100_w1)
  register_method("LRC-BART-f90-w0.9", fit_lrcbart_sa_f90_w09)
  register_method("LRC-BART-100-w0.9", fit_lrcbart_sa_100_w09)
  register_method("LRC-BART-full", fit_lrcbart_sa_full)
  invisible(TRUE)
}

single_arm_methods <- function(outcome) {
  c(if (outcome == "gaussian") c("LM-CP", "HierLM") else c("AFT-CP", "HierAFT"),
    "BART-CP", "LRC-BART-f90-w1", "LRC-BART-100-w1", "LRC-BART-f90-w0.9",
    "LRC-BART-100-w0.9", "LRC-BART-full")
}
