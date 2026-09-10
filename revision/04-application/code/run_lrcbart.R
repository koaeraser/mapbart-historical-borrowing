# run_lrcbart.R: LRC-BART in single-arm mode on the EloKRd + UCMM data
# (method_spec.md Sections A1, A3, C, F). f is AFT-BART on the 253 UCMM
# controls, g a prior draw (no RCT controls inform it), f_1 = f + g at the 30
# EloKRd profiles; the treated arm is plain AFT-BART on EloKRd (50 trees, as in
# the paper). Configurations: H_f in {10 (the paper's), 50}, w fixed at 1
# (reported) or 0.9 (sensitivity), H_g = 5 with split prior (0.5, 3) and
# nu0 = 3 (Phase 2b defaults), ESS_0 calibration on UCMM at N_target = 30 x m,
# m in {1, 0.75, 0.5}, and at the fractions 0.9 / 0.5 / 0.25 of the ceiling
# sigma_1^2 / V_mu^f. Calibration once per (outcome, coding, H_f) from a
# 2000 + 4000 Stage-1 chain; then 4 independent chains of 2000 warm-up + 2000
# draws for the control and treated fits. Estimands from lrc_ate (closed-form
# log-normal RMST at tau = 5 years, standardized to the 30 profiles; median
# ratio secondary). Convergence: split-Rhat and bulk/tail ESS (posterior
# package) on log(RMST ratio) across the 4 chains. One core.
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({library(lrcbart); library(posterior)})
APP <- file.path(Sys.getenv("LRC_REVISION_ROOT", "revision"), "04-application")
source(file.path(APP, "code/prep_data.R"))
TAU <- 5; H_G <- 5; N_BURN <- 2000; N_DRAW <- 2000; N_CHAINS <- 4
targets <- list(m1 = 30, m075 = 22.5, m05 = 15, f90 = c("frac", 0.9), f50 = c("frac", 0.5), f25 = c("frac", 0.25))
configs <- expand.grid(outcome = c("PFS", "OS"), coding = "harm", H_f = c(10, 50), w = c(1, 0.9), stringsAsFactors = FALSE)
configs <- rbind(configs, data.frame(outcome = c("PFS", "OS"), coding = "orig", H_f = 10, w = 1))

rows <- list(); draws_keep <- list(); maps <- list(); gsum <- list(); calib_keep <- list()
trt_cache <- list()
t_all <- Sys.time()
for (ci in seq_len(nrow(configs))) {
  cf <- configs[ci, ]
  dat <- load_mm(cf$outcome)
  X <- as.matrix(if (cf$coding == "orig") dat$X_orig else dat$X_harm)
  ti <- dat$trt_idx; ui <- dat$ctrl_idx
  x_rct <- X[ti, , drop = FALSE]; x_rwd <- X[ui, , drop = FALSE]
  y_rwd <- log(dat$time[ui]); st_rwd <- dat$status[ui]
  # treated arm, cached per (outcome, coding, chain)
  for (ch in seq_len(N_CHAINS)) {
    key <- paste(cf$outcome, cf$coding, ch)
    if (is.null(trt_cache[[key]]))
      trt_cache[[key]] <- bart_fit(NULL, x_rct, status = dat$status[ti], time = dat$time[ti],
                                   H_f = 50, n_burn = N_BURN, n_draw = N_DRAW, seed = 1000L + ch)
  }
  # Stage-1 fit and calibration, once per (outcome, coding, H_f)
  ckey <- paste(cf$outcome, cf$coding, cf$H_f)
  if (is.null(calib_keep[[ckey]])) {
    st1 <- bart_fit(NULL, x_rwd, status = st_rwd, time = dat$time[ui], H_f = cf$H_f,
                    n_burn = N_BURN, n_draw = 2 * N_DRAW, seed = 1L)
    cal_list <- lapply(targets, function(tg)
      lrc_ess_calibrate(y_rwd, x_rwd, x_rct, N_target = tg, status_rwd = st_rwd, fit_stage1 = st1,
                        H_g = H_G, nu0 = 3, H_f = cf$H_f, n_burn = N_BURN, n_draw = 2 * N_DRAW, seed = 1L, warn = FALSE))
    calib_keep[[ckey]] <- cal_list
  }
  cal_list <- calib_keep[[ckey]]
  for (tg in names(targets)) {
    cal <- cal_list[[tg]]
    t0 <- Sys.time()
    ate_r <- ate_m <- matrix(NA_real_, N_DRAW, N_CHAINS)
    g_all <- NULL; ess_real <- numeric(N_CHAINS); tau0 <- numeric(N_CHAINS); acc <- NULL
    for (ch in seq_len(N_CHAINS)) {
      ctrl <- lrc_bart_aft(dat$time[ui], st_rwd, x_rwd, rep(2L, length(ui)),
                           H_f = cf$H_f, H_g = H_G, alpha_g = 0.5, beta_g = 3, nu0 = 3,
                           s0_sq = cal$s0_sq, calib = cal, w = cf$w,
                           n_burn = N_BURN, n_draw = N_DRAW, seed = 100L * ch + ci, single_arm = TRUE)
      a <- lrc_ate(trt_cache[[paste(cf$outcome, cf$coding, ch)]], ctrl, x_rct, outcome = "survival", tau = TAU)
      ate_r[, ch] <- a$ate_rmst; ate_m[, ch] <- a$ate_med
      g_all <- rbind(g_all, lrc_discrepancy(ctrl, x_rct))
      ess_real[ch] <- lrc_ess_realized(ctrl, x_rct)$mean
      tau0[ch] <- mean(ctrl$tau0_sq)
      if (ch == 1) { acc <- ctrl$accept
        maps[[paste(cf$outcome, cf$coding, cf$H_f, cf$w, tg)]] <-
          lrc_ess_map(cal$fit_stage1, x_rct, ctrl$tau0_sq, H_g = H_G, sigma1_sq = cal$sigma1_sq) }
    }
    gsum[[paste(cf$outcome, cf$coding, cf$H_f, cf$w, tg)]] <-
      data.frame(profile = seq_len(nrow(x_rct)), mean = colMeans(g_all), sd = apply(g_all, 2, sd),
                 lo = apply(g_all, 2, quantile, 0.025), hi = apply(g_all, 2, quantile, 0.975))
    draws_keep[[paste(cf$outcome, cf$coding, cf$H_f, cf$w, tg)]] <- list(rmst = ate_r, med = ate_m)
    diag <- function(m) { lm <- log(m); c(rhat = rhat(lm), ess_bulk = ess_bulk(lm), ess_tail = ess_tail(lm)) }
    for (es in c("rmst", "med")) {
      m <- if (es == "rmst") ate_r else ate_m; v <- as.vector(m); d <- diag(m)
      rows[[length(rows) + 1]] <- data.frame(
        outcome = cf$outcome, coding = cf$coding, H_f = cf$H_f, H_g = H_G, w = cf$w, target = tg,
        N_target = cal$N_target, ceiling = cal$ceiling, ess_capped = cal$ess_capped, s0_sq = cal$s0_sq,
        ess0 = cal$ess0_at_s0, ess_realized = mean(ess_real), V_mu_f = cal$V_mu_f, sigma1_sq = cal$sigma1_sq,
        c_g = cal$c_g, tau0_sq_post = mean(tau0), g_sd_mean = mean(gsum[[length(gsum)]]$sd),
        estimand = es, mean = mean(v), median = median(v), lo = unname(quantile(v, 0.025)),
        hi = unname(quantile(v, 0.975)), p_gt1 = mean(v > 1),
        rhat = d["rhat"], ess_bulk = d["ess_bulk"], ess_tail = d["ess_tail"],
        acc_f_grow = acc["f", "grow"], acc_g_grow = acc["g", "grow"],
        secs = as.numeric(Sys.time() - t0, units = "secs"))
    }
    cat(sprintf("%s %s H_f=%d w=%g %s: N=%.1f ceil=%.1f s0=%.2g ess0=%.1f rmst %.3f (%.3f, %.3f) rhat %.3f  [%.0f s]\n",
                cf$outcome, cf$coding, cf$H_f, cf$w, tg, cal$N_target, cal$ceiling, cal$s0_sq, cal$ess0_at_s0,
                mean(ate_r), quantile(ate_r, 0.025), quantile(ate_r, 0.975), diag(ate_r)["rhat"],
                as.numeric(Sys.time() - t0, units = "secs")))
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  write.csv(tab, file.path(APP, "results/lrcbart_summary.csv"), row.names = FALSE)
  saveRDS(list(draws = draws_keep, maps = maps, gsum = gsum,
               calib = lapply(calib_keep, function(cl) lapply(cl, function(cal) cal[setdiff(names(cal), "fit_stage1")]))),
          file.path(APP, "results/lrcbart_draws.rds"))
}
cat("total", format(Sys.time() - t_all), "\n")
