# =============================================================================
# RMST(tau)-ratio helpers for the survival single-arm borrowing methods.
#
# The default estimand in these scripts is the population MEDIAN-survival ratio,
# computed from the mixture-of-lognormals survival curve S_mix(t) = mean_i
# Sbar((log t - mu_i)/sigma).  These helpers compute the analogous RESTRICTED
# MEAN survival time (RMST) ratio, RMST(tau) = integral_0^tau S_mix(t) dt, as an
# ADDITIONAL output.  RMST is a bounded, well-conditioned functional (no root
# finding, no 30-cap / degenerate-draw clamping needed) but its true value is
# tau-dependent and varies by replicate, so bias is measured against the
# per-replicate true RMST(tau)-ratio (computed below from the DGP truth).
#
# All inputs use the SAME posterior objects the median computation uses:
#   mu_*  : n_draws x n_test matrix of posterior log-survival means (yhat/f.test)
#   sig_* : length-n_draws vector of posterior residual SDs
# =============================================================================

# population mixture-of-lognormals MEDIAN survival, one value per posterior draw.
# Mirrors the inline root-finding in the method scripts (lower floor 0.05, cap 30,
# degenerate draws clamped to the cap). `sig_vec` is recycled to nrow(mu_mat) so the
# control surface draws can be paired with the TREATMENT model's sigma draws even
# when the two chains have different lengths.
pop_median_draws <- function(mu_mat, sig_vec, lower = 0.05, cap = 30) {
  w <- rep(1 / ncol(mu_mat), ncol(mu_mat))
  Smix <- function(t, mu, sig) sum(w * pnorm(log(t), mu, sig, lower.tail = FALSE))
  n <- nrow(mu_mat)
  sig_vec <- rep_len(sig_vec, n)
  upper <- pmin(exp(apply(mu_mat, 1, max) + 6 * sig_vec), cap)
  out <- numeric(n)
  for (m in seq_len(n)) {
    f <- function(t) Smix(t, mu_mat[m, ], sig_vec[m]) - 0.5
    if (f(upper[m]) > 0) { out[m] <- upper[m]; next }      # degenerate -> clamp to cap
    lo <- max(1e-10, lower)
    out[m] <- if (f(lo) < 0) lo else uniroot(f, c(lo, upper[m]))$root
  }
  out
}

# population RMST(tau) for one posterior draw: integrate the mixture survival
# curve over [0, tau] via a Riemann sum on a fixed grid.
.rmst_pop_grid <- function(mu, sig, lg, tau) {
  # mu: length-n log-survival means; sig: scalar; lg: log(grid) of length K
  mean(rowMeans(pnorm(outer(lg, mu, "-") / sig, lower.tail = FALSE))) * tau
}

# per-draw population RMST(tau) for ONE arm (length-n_draws vector)
rmst_draws <- function(mu_mat, sig_vec, tau = 3, K = 120) {
  nd <- min(nrow(mu_mat), length(sig_vec))
  lg <- log(seq(1e-3, tau, length.out = K))
  vapply(seq_len(nd), function(m) .rmst_pop_grid(mu_mat[m, ], sig_vec[m], lg, tau), numeric(1))
}

# per-(draw, subject) RMST(tau): n_draws x n_subjects matrix. The population RMST
# per draw is rowMeans() of this (mean over subjects), so one call gives both the
# population and the subject-level RMST at no extra cost.
rmst_subj_mat <- function(mu_mat, sig_vec, tau = 3, K = 120) {
  nd <- min(nrow(mu_mat), length(sig_vec))
  lg <- log(seq(1e-3, tau, length.out = K))
  t(vapply(seq_len(nd), function(m)
    colMeans(pnorm(outer(lg, mu_mat[m, ], "-") / sig_vec[m], lower.tail = FALSE)) * tau,
    numeric(ncol(mu_mat))))
}

# per-draw posterior RMST(tau)-ratio (treatment / control).
# `floor` bounds the control RMST away from 0: at weak borrowing a few posterior
# draws give a near-zero control survival (RMST_ctrl ~ 0), which makes the ratio
# (and its var-based SD) blow up. This mirrors the median estimand's `lower=0.05`
# floor on the control median; the quantile CI is unaffected, only the SD is.
rmst_ratio_post <- function(mu_t, sig_t, mu_c, sig_c, tau = 3, K = 120, floor = 0.05) {
  rt <- rmst_draws(mu_t, sig_t, tau, K); rc <- rmst_draws(mu_c, sig_c, tau, K)
  n <- min(length(rt), length(rc))
  rt[seq_len(n)] / pmax(rc[seq_len(n)], floor)
}

# true (DGP) population RMST(tau) per arm for one replicate: c(trt=, ctrl=).
# Uses the values saved by data_gen (true_rmst_{trt,ctrl}_pop) when present and the
# restriction horizon matches (true_rmst_tau == tau); otherwise computes on the fly
# from lp/eff_i/sigma_rct (same mixture the DGP uses for the median).
true_rmst_arms <- function(data_tmp, tau = 3, K = 120) {
  if (!is.null(data_tmp$true_rmst_trt_pop) && !is.null(data_tmp$true_rmst_tau) &&
      isTRUE(all.equal(data_tmp$true_rmst_tau, tau))) {
    return(c(trt = unname(data_tmp$true_rmst_trt_pop), ctrl = unname(data_tmp$true_rmst_ctrl_pop)))
  }
  rct <- data_tmp$X$D == 1
  MUt <- data_tmp$lp[rct]
  MUc <- MUt - data_tmp$eff_i
  sig <- data_tmp$sigma_rct
  lg  <- log(seq(1e-3, tau, length.out = K))
  c(trt = .rmst_pop_grid(MUt, sig, lg, tau), ctrl = .rmst_pop_grid(MUc, sig, lg, tau))
}

# true (DGP) PER-SUBJECT RMST(tau)-ratio for the D==1 population (length n_test).
# RMST_i(mu) = integral_0^tau Sbar((log t - mu_i)/sigma) dt for each subject i.
true_rmst_subject <- function(data_tmp, tau = 3, K = 120) {
  MUt <- log(data_tmp$true_median_trt)   # per-subject treated log-survival mean
  MUc <- log(data_tmp$true_median_ctrl)  # per-subject control  log-survival mean
  sig <- data_tmp$sigma_rct
  lg  <- log(seq(1e-3, tau, length.out = K))
  rmst1 <- function(mu) colMeans(pnorm(outer(lg, mu, "-") / sig, lower.tail = FALSE)) * tau
  rmst1(MUt) / rmst1(MUc)
}

# true (DGP) population RMST(tau)-ratio for the D==1 population of one replicate.
# effect = NULL  -> the actual treatment effect (uses saved/true arm values)  [true value]
# effect = scalar -> counterfactual where treated = control + effect          [eff_star]
true_rmst_ratio <- function(data_tmp, tau = 3, K = 120, effect = NULL) {
  if (is.null(effect)) {
    a <- true_rmst_arms(data_tmp, tau, K)
    return(unname(a["trt"] / a["ctrl"]))
  }
  rct <- data_tmp$X$D == 1
  MUc <- data_tmp$lp[rct] - data_tmp$eff_i              # control counterfactual (per subject)
  MUt <- MUc + effect
  sig <- data_tmp$sigma_rct
  lg  <- log(seq(1e-3, tau, length.out = K))
  .rmst_pop_grid(MUt, sig, lg, tau) / .rmst_pop_grid(MUc, sig, lg, tau)
}

# one-row list of RMST(tau)-ratio metrics, mirroring the median-ratio metrics.
# control_sigma selects the sigma used for the CONTROL-arm RMST (the treatment arm
# always uses its own sig_t, which is fine since its RMST saturates near tau):
#   "trt"      sig_t (treatment model's sigma-hat). Best at large n; at small n the
#              trt sigma-hat is inflated (e.g. 1.9 vs true 0.6 at n=30), which
#              over-inflates the control RMST and biases the ratio DOWN.
#   "own"      sig_c (control model's own sigma-hat; RWD, large-sample, stable).
#   "adaptive" pick whichever arm has the smaller central sigma-hat (the less-inflated
#              one): trt at large n, control at small n. Robust across sample sizes.
# NOTE: pass the CONTROL model's own sigma-hat as sig_c for "own"/"adaptive" to work.
compute_rmst_metrics <- function(mu_t, sig_t, mu_c, sig_c, data_tmp,
                                 tau = 3, K = 120, threshold = 0.95, floor = 0.05,
                                 control_sigma = "trt") {
  sig_ctrl <- switch(control_sigma,
                     own      = sig_c,
                     adaptive = if (median(sig_t, na.rm = TRUE) <= median(sig_c, na.rm = TRUE)) sig_t else sig_c,
                     sig_t)   # "trt" / default
  true_sig_c <- if (!is.null(data_tmp$sigma_rwd)) data_tmp$sigma_rwd else data_tmp$sigma_rct
  Rt <- rmst_subj_mat(mu_t, sig_t,    tau, K)  # treatment arm: always its own sigma-hat
  Rc <- rmst_subj_mat(mu_c, sig_ctrl, tau, K)  # control arm: per control_sigma
  n  <- min(nrow(Rt), nrow(Rc)); Rt <- Rt[seq_len(n), , drop = FALSE]; Rc <- Rc[seq_len(n), , drop = FALSE]
  rt <- rowMeans(Rt); rc <- rowMeans(Rc)       # population RMST per draw (mean over subjects)
  ps <- rt / pmax(rc, floor)                   # population ratio (control RMST floored away from 0)

  eff      <- true_rmst_ratio(data_tmp, tau, K)                                  # per-rep ratio truth
  eff_star <- true_rmst_ratio(data_tmp, tau, K, effect = data_tmp$treat_eff_star)
  arms     <- true_rmst_arms(data_tmp, tau, K)                                   # per-arm truth

  delta <- median(ps, na.rm = TRUE)
  qs    <- quantile(ps, c(0.025, 0.975), na.rm = TRUE)
  dec   <- mean(ps > eff_star, na.rm = TRUE)

  # subject-specific RMST(tau)-ratio (per-subject estimate vs per-subject truth)
  delta_i <- apply(Rt / pmax(Rc, floor), 2, median, na.rm = TRUE)   # n_subjects
  eff_i   <- true_rmst_subject(data_tmp, tau, K)                    # per-subject truth

  list(rmst_tau      = tau,
       rmst_true     = eff,
       rmst_hat      = delta,
       bias_rmst     = delta - eff,
       sd_rmst       = sqrt(var(ps, na.rm = TRUE)),
       ci_rmst       = unname(diff(qs)),
       coverage_rmst = unname(qs[1] <= eff & qs[2] >= eff),
       # additional within/across-replicate metrics, mirroring the median-ratio ones
       rmse_rmst       = (delta - eff)^2,                          # across-rep: sqrt(mean(.)) in plot
       w1distance_rmst = mean(abs(ps - eff), na.rm = TRUE),        # E|post - truth|
       w2distance_rmst = sqrt(mean((ps - eff)^2, na.rm = TRUE)),   # sqrt E(post - truth)^2
       decision_rmst = dec,
       tp_rmst            = as.numeric(dec > 0.95),       # power (uncalibrated, mirrors median tp)
       tp_calibrated_rmst = as.numeric(dec > threshold),  # calibrated power (mirrors median tp_calibrated)
       # arm-specific population RMST: estimate (median over draws) vs DGP truth,
       # mirroring the *.median.pop metrics (uses UNFLOORED RMST draws).
       "bias.trt.rmst.pop"        = median(rt, na.rm = TRUE) - unname(arms["trt"]),
       "sd.trt.rmst.pop"          = sqrt(var(rt, na.rm = TRUE)),   # posterior SD (matches ATE sd_rmst)
       "w2distance.trt.rmst.pop"  = sqrt(mean((rt - unname(arms["trt"]))^2, na.rm = TRUE)),
       "bias.ctrl.rmst.pop"       = median(rc, na.rm = TRUE) - unname(arms["ctrl"]),
       "sd.ctrl.rmst.pop"         = sqrt(var(rc, na.rm = TRUE)),   # posterior SD (matches ATE sd_rmst)
       "w2distance.ctrl.rmst.pop" = sqrt(mean((rc - unname(arms["ctrl"]))^2, na.rm = TRUE)),
       # subject-wise RMST: bias and PEHE over subjects (mirrors bias_subj / pehe_subj)
       bias_subj_rmst = mean(delta_i - eff_i, na.rm = TRUE),
       pehe_subj_rmst = sqrt(mean((delta_i - eff_i)^2, na.rm = TRUE)),
       # control-arm residual-SD (sigma) metrics, mirroring the treatment-arm
       # bias.trt.sigma / sd.trt.sigma / w2distance.trt.sigma (uses the control
       # model's OWN sigma-hat draws sig_c, vs the true control sigma).
       "bias.ctrl.sigma"       = mean(sig_c, na.rm = TRUE) - true_sig_c,
       "sd.ctrl.sigma"         = sd(sig_c, na.rm = TRUE),
       "w2distance.ctrl.sigma" = sqrt(mean((sig_c - true_sig_c)^2, na.rm = TRUE)))
}
