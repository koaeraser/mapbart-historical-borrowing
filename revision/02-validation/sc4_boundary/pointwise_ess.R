# pointwise_ess.R: pointwise prior and realized ESS at covariate profiles for
# LRC-BART, built from the package's own Stage-1 / Stage-2 machinery without
# touching the sampler.
#
# Estimand-level definitions (method_spec.md Section C; lrcbart/R/lrcbart.R):
#   ESS_0        = sigma_1^2 E_{tau0^2 ~ Inv-chi2(nu0, s0^2)}[1 / (V_mu^f + c_g H_g tau0^2)]
#   realized ESS = E_post[sigma_1^2 / (V_mu^f + V_g)], V_g the estimand-level
#                  discrepancy variance given the posterior (partition, z, tau^2)
# with mu = N^-1 sum_i f_1(x_i) the RCT-standardized control mean, V_mu^f its
# Stage-1 (RWD-only BART) posterior variance, and c_g the prior probability
# that two RCT profiles share a g leaf.
#
# Pointwise versions at a profile x (spec Section C, last paragraph, with the
# expectation over the calibrated prior of tau0^2 instead of a plug-in):
#   ESS_0(x)        = sigma_1^2 E_{tau0^2}[1 / (V_f(x) + H_g tau0^2)]
#   ceiling(x)      = sigma_1^2 / V_f(x)                       (tau0^2 -> 0)
#   realized ESS(x) = E_post[sigma_1^2 / (V_f(x) + sum_h tau^2_{z_h(x)})]
#   post ESS(x)     = sigma_1^2 / Var_post(f_1(x))             (literal version)
# where V_f(x) = Var(f(x) | RWD) from the whole Stage-1 chain, sigma_1^2 is
# the calibration's RCT residual variance (prior versions) or the posterior
# draws of sigma_1^2 (realized version), and z_h(x) is the spike/slab state of
# the g leaf of tree h that contains x. The pointwise leaf-sharing constant is
# 1 (x shares its own leaf), which is why c_g does not appear: every one of
# the H_g trees contributes a full tau^2 to the variance at a single point,
# whereas at the estimand level only the fraction c_g of profile pairs share
# a leaf. Because Var of a mean is below the mean of the variances (positive
# but imperfect correlation of f_1 across profiles), the pointwise ESS is
# smaller than ESS_0 at every profile; the two are consistent, not equal.

suppressPackageStartupMessages(library(lrcbart))

# Refit the calibration and the control model exactly as fit_lrcbart() does
# (01-code/R/fit_lrcbart.R, steps 1 and 2, warm_start = FALSE), consuming the
# RNG in the same order, so that with set.seed(fit_seed(seed, r, method))
# before the call the fit is the one of the full study. The treated-arm fit
# (step 3) is skipped: it does not enter the control surface or the ESS.
lrc_refit_parts <- function(dat, N_target = 100, w_prior = c(1, 1), w_fixed = NULL,
                            H_f = 50, H_g = 5, alpha_g = 0.5, beta_g = 3, nu0 = 3, q = 0.95,
                            n_burn = 1000L, n_draw = 1000L, calib_grid = NULL,
                            tau1_mult = 1, n_min_g = 5, g_sweeps = 1) {
  xc <- xmat(dat$rct_ctrl); xw <- xmat(dat$rwd_ctrl)
  x_rct <- dat$x_rct
  src <- rep(c(1L, 2L), c(nrow(xc), nrow(xw)))
  y_rwd <- dat$rwd_ctrl$y
  seed <- sample.int(.Machine$integer.max, 1)  # identical to fit_lrcbart()
  cal <- lrc_ess_calibrate(y_rwd, xw, x_rct, N_target = N_target, grid = calib_grid, q = q,
                           status_rwd = NULL,
                           y_rct_ctrl = dat$rct_ctrl$y, x_rct_ctrl = xc, status_rct_ctrl = NULL,
                           nu0 = nu0, H_g = H_g, alpha_g = alpha_g,
                           beta_g = beta_g, H_f = H_f, n_burn = n_burn, n_draw = n_draw,
                           seed = seed, warn = FALSE)
  ctrl <- lrc_bart(c(dat$rct_ctrl$y, dat$rwd_ctrl$y), rbind(xc, xw), src,
                   H_f = H_f, H_g = H_g, alpha_g = alpha_g, beta_g = beta_g, nu0 = nu0,
                   s0_sq = cal$s0_sq, calib = cal, w = w_fixed, a_w = w_prior[1], b_w = w_prior[2],
                   tau1_mult = tau1_mult, n_min_g = n_min_g, g_sweeps = g_sweeps, g_init = NULL,
                   n_burn = n_burn, n_draw = n_draw, seed = seed + 1L)
  list(cal = cal, ctrl = ctrl, xc = xc, seed = seed)
}

# Draws of tau0^2 from the calibrated scaled-Inv-chi2(nu0, s0^2) prior.
tau0_prior_draws <- function(cal, n_tau = 4000, seed = 1) {
  set.seed(seed)
  cal$nu0 * cal$s0_sq / stats::rchisq(n_tau, cal$nu0)
}

# Pointwise prior ESS at the profiles x (rows). Returns a data.frame with
# V_f(x), ESS_0(x), ceiling(x). Vectorized over the tau0^2 draws.
pointwise_prior_ess <- function(cal, x, tau0 = tau0_prior_draws(cal)) {
  f <- predict(cal$fit_stage1, x)               # n_draw x n_x Stage-1 draws
  V_f <- apply(f, 2, stats::var)                # whole-chain variance, as V_mu^f
  H_g <- cal$H_g; s1 <- cal$sigma1_sq
  ess0 <- vapply(V_f, function(v) s1 * mean(1 / (v + H_g * tau0)), 0)
  data.frame(V_f = V_f, ess0_x = ess0, ceiling_x = s1 / V_f)
}

# Estimand-level ESS_0 recomputed with the whole-chain V_mu^f (the package's
# reported ess0_at_s0 averages sigma_1^2 / (V_b + ...) over the 20 rescaled
# blocks, which is larger by Jensen).
ess0_whole_chain <- function(cal, tau0 = tau0_prior_draws(cal)) {
  cal$sigma1_sq * mean(1 / (cal$V_mu_f + cal$c_g * cal$H_g * tau0))
}

# Pointwise realized ESS and the literal posterior version at the profiles x.
pointwise_realized_ess <- function(ctrl, cal, x, V_f) {
  pf <- lrcbart:::predict_forest_cpp(ctrl$draws$g, lrcbart:::as_x_matrix(x))
  ns <- pf$nspike                                # n_draw x n_x spike trees at x
  Vg <- ns * ctrl$tau0_sq + (ctrl$H_g - ns) * ctrl$tau1_sq   # recycled by column
  ess_real <- colMeans(ctrl$sigma1_sq / sweep(Vg, 2, V_f, "+"))
  f1 <- predict(ctrl, x, arm = "rct_control")
  ess_post <- mean(ctrl$sigma1_sq) / apply(f1, 2, stats::var)
  data.frame(ess_real_x = ess_real, ess_post_x = ess_post, nspike_x = colMeans(ns),
             surface_x = colMeans(f1), g_mean_x = colMeans(pf$value))
}

# Regional ESS: the paper's estimand-level definition applied to the
# region-standardized control mean mu_S = |S|^-1 sum_{i in S} f_1(x_i) over the
# RCT profiles x_S in a region S (R or R^c):
#   ESS_0(S)        = sigma_1^2 E_{tau0^2}[1 / (V_{mu_S}^f + c_g(S) H_g tau0^2)]
#   realized ESS(S) = E_post[sigma_1^2 / (V_{mu_S}^f + V_g(S))]
# with V_{mu_S}^f the Stage-1 variance of mu_S, c_g(S) the leaf-sharing
# constant on x_S under the g-tree prior (lrc_cg) and V_g(S) the exact
# estimand-level discrepancy variance on x_S given the posterior partition
# (forest_gvar_cpp, as in lrc_ess_realized). With S = all RCT profiles this
# is ESS_0 / the realized ESS of the package (whole-chain V).
regional_ess <- function(cal, ctrl, x_S, x_grid, tau0 = tau0_prior_draws(cal), n_cg = 2000, seed_cg = 1) {
  f <- predict(cal$fit_stage1, x_S)
  V_S <- stats::var(rowMeans(f))
  c_gS <- lrc_cg(x_S, x_grid = x_grid, alpha_g = ctrl$hyper$alpha_g,
                 beta_g = ctrl$hyper$beta_g, n_sim = n_cg, seed = seed_cg)
  ess0 <- cal$sigma1_sq * mean(1 / (V_S + c_gS * cal$H_g * tau0))
  Vg <- lrcbart:::forest_gvar_cpp(ctrl$draws$g, lrcbart:::as_x_matrix(x_S), ctrl$tau0_sq, ctrl$tau1_sq)
  real <- mean(ctrl$sigma1_sq / (V_S + Vg))
  data.frame(n_S = nrow(x_S), V_S = V_S, c_g_S = c_gS, ess0_S = ess0, ceiling_S = cal$sigma1_sq / V_S, ess_real_S = real)
}
