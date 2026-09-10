# fit_lrcbart.R: harness wrapper for the proposed method (package lrcbart,
# revision/01-code/lrcbart). Same interface and return object as the
# comparators in comparators.R: fit_<name>(dat, outcome, ...) returning
# make_fit()'s list (draws, surface, map, ess, extra).
#
# Pipeline per replicate (method_spec.md v2, Sections B and C):
#   1. Stage-1 standard BART on the RWD controls alone, evaluated at all
#      n_rct RCT profiles; ESS_0 calibration of s0^2 by the Pr-rule at
#      N_target (q = 0.95, nu0 = 3, c_g by Monte Carlo under the g-tree prior).
#   2. LRC-BART on the pooled controls (RCT controls source 1, RWD source 2),
#      H_f = 50, H_g = 5, split prior (0.5, 3), w ~ Beta(1, 1) unless fixed,
#      tau1^2 by the range rule, tau0^2 learned from the spike leaves (the
#      Phase 2b configuration, lrcbart/NOTES.md "Phase 2b tuning"; H_g = 10
#      with w ~ Beta(4, 1) was the development-run configuration).
#   3. Plain BART on the RCT treated subjects (same f code path, one source).
#   4. Standardized estimands over the RCT profiles through the harness's
#      standardized_draws(), so the estimand code is shared with the comparators.
# The compatibility map B(x) = P(|g(x)| < c | data), c = map_c (0.5 on the
# outcome scale for Gaussian, on the log-time scale for survival), is reported
# at the RCT control profiles (the harness's map_R, map_Rc, map_all); the
# ensemble-level all-spike map is kept in extra$map_spike. ess = list(prior =
# ESS_0 at the calibrated s0^2, realized = posterior realized ESS, ceiling,
# ess_capped); the calibration's sigma_1^2 is the RCT-control residual variance.
#
# Requires: library(lrcbart) installed from revision/01-code/lrcbart
# (R CMD INSTALL or devtools::install). Sourced after comparators.R so that
# register_method() is available; if it is not, the functions are defined
# without registration.

suppressPackageStartupMessages(library(lrcbart))

# Tuning arguments (Phase 2b): tau1_mult multiplies the range-rule slab scale,
# n_min_g is the g leaf minimum size, g_sweeps the number of g-ensemble sweeps
# per iteration, warm_start = TRUE initializes tree 1 of g from a greedy tree
# (depth 2) fitted to the BART-PP-implied source offset at the RCT control
# profiles (diagnostic for mixing). N_target may be an absolute ESS or
# c("frac", p) for the fraction p of the ESS ceiling.
fit_lrcbart <- function(dat, outcome = c("gaussian", "survival"), N_target = 100,
                        w_prior = c(1, 1), w_fixed = NULL, H_f = 50, H_g = 5,
                        alpha_g = 0.5, beta_g = 3, nu0 = 3, q = 0.95,
                        n_burn = N_BURN, n_draw = N_DRAW, calib_grid = NULL,
                        map_c = 0.5, tau1_mult = 1, n_min_g = 5, g_sweeps = 1,
                        warm_start = FALSE, ...) {
  outcome <- match.arg(outcome)
  t0 <- Sys.time()
  xt <- xmat(dat$rct_trt); xc <- xmat(dat$rct_ctrl); xw <- xmat(dat$rwd_ctrl)
  x_rct <- dat$x_rct
  src <- rep(c(1L, 2L), c(nrow(xc), nrow(xw)))
  gauss <- outcome == "gaussian"
  y_rwd <- if (gauss) dat$rwd_ctrl$y else log(dat$rwd_ctrl$time)
  status_rwd <- if (gauss) NULL else dat$rwd_ctrl$status
  seed <- sample.int(.Machine$integer.max, 1)  # from the harness's set.seed()
  # 0. optional warm start: BART with the source indicator as a covariate
  #    (the BART-PP model) on the pooled controls, offset D = 1 minus D = 0 at
  #    the RCT control profiles, then a greedy depth-2 tree on that offset
  g_init <- NULL
  if (isTRUE(warm_start) && nrow(xc) > 0) {
    xd <- cbind(rbind(xc, xw), D = as.numeric(src == 1L))
    yd <- if (gauss) c(dat$rct_ctrl$y, dat$rwd_ctrl$y) else log(c(dat$rct_ctrl$time, dat$rwd_ctrl$time))
    bpp <- if (gauss) bart_fit(yd, xd, H_f = H_f, n_burn = n_burn, n_draw = n_draw, seed = seed + 3L)
           else bart_fit(NULL, xd, status = c(dat$rct_ctrl$status, dat$rwd_ctrl$status),
                         time = exp(yd), H_f = H_f, n_burn = n_burn, n_draw = n_draw, seed = seed + 3L)
    off <- colMeans(predict(bpp, cbind(xc, D = 1))) - colMeans(predict(bpp, cbind(xc, D = 0)))
    g_init <- lrc_g_init_greedy(xc, off, H_g = H_g, maxdepth = 2, n_min = max(n_min_g, 2), c_slab = map_c)
  }

  # 1. calibration; sigma_1^2 is the posterior mean residual variance of a
  #    plain BART fitted to the RCT controls whenever they exist (see
  #    lrc_ess_calibrate, sigma1_method = "bart"), so the prior ESS_0 and the
  #    realized ESS are in the same (RCT) units
  y_rct_ctrl <- if (nrow(xc) == 0) NULL else if (gauss) dat$rct_ctrl$y else log(dat$rct_ctrl$time)
  status_rct_ctrl <- if (nrow(xc) == 0 || gauss) NULL else dat$rct_ctrl$status
  cal <- lrc_ess_calibrate(y_rwd, xw, x_rct, N_target = N_target, grid = calib_grid, q = q,
                           status_rwd = status_rwd,
                           y_rct_ctrl = y_rct_ctrl, x_rct_ctrl = xc, status_rct_ctrl = status_rct_ctrl,
                           nu0 = nu0, H_g = H_g, alpha_g = alpha_g,
                           beta_g = beta_g, H_f = H_f, n_burn = n_burn, n_draw = n_draw,
                           seed = seed, warn = FALSE)  # the cap is recorded in extra$ess_capped
  # 2. controls
  ctrl <- if (gauss) {
    lrc_bart(c(dat$rct_ctrl$y, dat$rwd_ctrl$y), rbind(xc, xw), src,
             H_f = H_f, H_g = H_g, alpha_g = alpha_g, beta_g = beta_g, nu0 = nu0,
             s0_sq = cal$s0_sq, calib = cal, w = w_fixed, a_w = w_prior[1], b_w = w_prior[2],
             tau1_mult = tau1_mult, n_min_g = n_min_g, g_sweeps = g_sweeps, g_init = g_init,
             n_burn = n_burn, n_draw = n_draw, seed = seed + 1L, ...)
  } else {
    lrc_bart_aft(c(dat$rct_ctrl$time, dat$rwd_ctrl$time), c(dat$rct_ctrl$status, dat$rwd_ctrl$status),
                 rbind(xc, xw), src,
                 H_f = H_f, H_g = H_g, alpha_g = alpha_g, beta_g = beta_g, nu0 = nu0,
                 s0_sq = cal$s0_sq, calib = cal, w = w_fixed, a_w = w_prior[1], b_w = w_prior[2],
                 tau1_mult = tau1_mult, n_min_g = n_min_g, g_sweeps = g_sweeps, g_init = g_init,
                 n_burn = n_burn, n_draw = n_draw, seed = seed + 1L, ...)
  }
  # 3. treated arm
  trt <- if (gauss) {
    bart_fit(dat$rct_trt$y, xt, H_f = H_f, n_burn = n_burn, n_draw = n_draw, seed = seed + 2L)
  } else {
    bart_fit(NULL, xt, status = dat$rct_trt$status, time = dat$rct_trt$time, H_f = H_f,
             n_burn = n_burn, n_draw = n_draw, seed = seed + 2L)
  }
  # 4. estimands at the RCT profiles
  mu1 <- predict(trt, x_rct)
  mu0 <- predict(ctrl, x_rct, arm = "rct_control")
  sig1 <- if (gauss) NULL else sqrt(trt$sigma1_sq)
  sig0 <- if (gauss) NULL else sqrt(ctrl$sigma1_sq)
  # compatibility map P(|g(x)| < c | data) at the RCT control profiles (the
  # harness's map_R / map_Rc); the ensemble-level all-spike map is kept in extra
  map <- lrc_borrowing_map(ctrl, xc, type = "abs_g", c = map_c)$mean
  map_spike <- lrc_borrowing_map(ctrl, xc, type = "all_spike")$mean
  ess_real <- lrc_ess_realized(ctrl, x_rct)
  gx <- lrc_discrepancy(ctrl, xc)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  make_fit(mu1, mu0, outcome, sig1, sig0, dat$constants$tau,
           map = map,
           ess = list(prior = cal$ess0_at_s0, realized = ess_real$mean,
                      ceiling = cal$ceiling, s0_sq = cal$s0_sq, c_g = cal$c_g,
                      V_mu_f = cal$V_mu_f, at_boundary = cal$at_boundary,
                      sigma1_sq_calib = cal$sigma1_sq, N_target = cal$N_target,
                      N_target_spec = cal$N_target_spec,
                      ess_capped = cal$ess_capped),
           extra = list(map_spike = map_spike, map_c = map_c, ess_capped = cal$ess_capped,
                        g_mean = colMeans(gx), g_sd = apply(gx, 2, sd),
                        config = list(H_g = H_g, alpha_g = alpha_g, beta_g = beta_g,
                                      tau1_mult = tau1_mult, n_min_g = n_min_g,
                                      g_sweeps = g_sweeps, warm_start = warm_start),
                        tau0_sq = mean(ctrl$tau0_sq), tau1_sq = mean(ctrl$tau1_sq),
                        w = mean(ctrl$w), spike_frac = mean(ctrl$K0 / pmax(ctrl$L_g, 1)),
                        L_g = mean(ctrl$L_g), accept = ctrl$accept,
                        sigma1 = sqrt(mean(ctrl$sigma1_sq)), sigma2 = sqrt(mean(ctrl$sigma2_sq)),
                        time = secs))
}

# Named variants matching the spec's comparator list (Section E). Absolute
# ESS targets 100 / 75 / 50 (capped at the ceiling sigma1^2 / V_mu^f when
# infeasible) and fractional targets f90 / f50 / f25 (90, 50 and 25 percent
# of the ceiling, always feasible); see Section A3 of the spec for the ladder.
fit_lrcbart_100 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 100, ...)
fit_lrcbart_75 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 75, ...)
fit_lrcbart_50 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 50, ...)
fit_lrcbart_f90 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = c("frac", 0.9), ...)
fit_lrcbart_f50 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = c("frac", 0.5), ...)
fit_lrcbart_f25 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = c("frac", 0.25), ...)
# sensitivity: the development-run g configuration (H_g = 10, w ~ Beta(4, 1))
fit_lrcbart_hg10 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 100, H_g = 10, w_prior = c(4, 1), ...)
# C-BART: w = 1 (all spike), same calibrated s0^2.
fit_cbart <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 100, w_fixed = 1, ...)
# free nonparametric offset: w = 0 (all slab), the tree analogue of BART-PP.
fit_lrcbart_w0 <- function(dat, outcome, ...) fit_lrcbart(dat, outcome, N_target = 100, w_fixed = 0, ...)

if (exists("register_method")) {
  register_method("LRC-BART-100", fit_lrcbart_100)
  register_method("LRC-BART-75", fit_lrcbart_75)
  register_method("LRC-BART-50", fit_lrcbart_50)
  register_method("LRC-BART-f90", fit_lrcbart_f90)
  register_method("LRC-BART-f50", fit_lrcbart_f50)
  register_method("LRC-BART-f25", fit_lrcbart_f25)
  register_method("LRC-BART-Hg10", fit_lrcbart_hg10)
  register_method("C-BART", fit_cbart)
  register_method("LRC-BART-w0", fit_lrcbart_w0)
}
