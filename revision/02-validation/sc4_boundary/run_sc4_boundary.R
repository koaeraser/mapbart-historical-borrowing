# run_sc4_boundary.R: Sc4 boundary metrics for the grader's first ask.
#
# (A) Control-surface RMSE at the RCT control profiles, overall, in R = {X5 > 2},
#     in R^c, and restricted to profiles more than 0.5 from the boundary
#     (R far: X5 > 2.5; R^c far: X5 < 1.5), for LRC-BART-100, BART-PP, BART-CP,
#     BART-NP. The full-study RDS files (02-validation/full/gaussian/results)
#     keep the posterior-mean surface at all 300 RCT profiles, so (A) is
#     computed from them without refitting; the replicate data are regenerated
#     with the study's seeds (rep_seed(2026, r)) for X5, the arm and the truth.
#     The all-profile version (the paper's Table 2 definition) is kept as a
#     cross-check column.
# (B) Pointwise prior ESS_0(x) (and the realized / posterior versions) at the
#     RCT control profiles for LRC-BART-100, from refits of the calibration and
#     the control model with the study's fit seeds (fit_seed(2026, r, method)),
#     see pointwise_ess.R. The refit surface is checked against the stored one.
#
# Scenarios: Sc1 (reference), Sc4_X5_d0.5, Sc4_X5_d1, Sc4_X5_d2, Sc5_d2;
# R = 200 (the first 200 of the study's 500; the (A) metrics are also written
# at R = 500 in summary_R500.csv / paired_R500.csv). OMP_NUM_THREADS = 1,
# mclapply with 12 cores. Usage: Rscript run_sc4_boundary.R [N_REP] [CORES].
Sys.setenv(OMP_NUM_THREADS = "1")
REV <- Sys.getenv("LRC_REVISION_ROOT", "revision")
CODE_ROOT <- file.path(REV, "01-code")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
HERE <- file.path(REV, "02-validation", "sc4_boundary")
source(file.path(HERE, "pointwise_ess.R"))
FULL <- file.path(REV, "02-validation", "full", "gaussian", "results")
OUT <- file.path(HERE, "results"); dir.create(OUT, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
N_REP <- if (length(args) >= 1) as.integer(args[1]) else 200L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 12L
N_ALL <- 500L
SEED <- 2026; METHOD <- "LRC-BART-100"
METHODS <- c("LRC-BART-100", "BART-PP", "BART-CP", "BART-NP")
SCEN <- c("Sc1", "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d2")
FAR <- 0.5
t_start <- Sys.time()

region_of <- function(x5) {
  data.frame(R = x5 > 2, Rc = x5 <= 2, R_far = x5 > 2 + FAR, Rc_far = x5 < 2 - FAR, all = TRUE)
}
rmse_by_region <- function(err, reg) {
  vapply(names(reg), function(k) if (any(reg[[k]])) sqrt(mean(err[reg[[k]]]^2)) else NA_real_, 0)
}

# ---------------------------------------------------------------------------
# (A) surface RMSE from the stored surfaces
# ---------------------------------------------------------------------------
surface_one <- function(sc, r) {
  dat <- gen_data_id(sc, "gaussian", seed = rep_seed(SEED, r), n_mc = 0)
  ctrl <- dat$z_rct == 0
  x5 <- dat$x_rct[, 5]
  reg_c <- region_of(x5[ctrl]); reg_all <- region_of(x5)
  rows <- lapply(METHODS, function(m) {
    f <- result_path(FULL, sc, m, r)
    if (!file.exists(f)) return(NULL)
    s <- readRDS(f)$surface
    err <- s - dat$mu0_rct
    rc <- rmse_by_region(err[ctrl], reg_c); ra <- rmse_by_region(err, reg_all)
    data.frame(scenario = sc, method = m, rep = r,
               rmse_all = rc[["all"]], rmse_R = rc[["R"]], rmse_Rc = rc[["Rc"]],
               rmse_R_far = rc[["R_far"]], rmse_Rc_far = rc[["Rc_far"]],
               rmse_allprof_all = ra[["all"]], rmse_allprof_R = ra[["R"]], rmse_allprof_Rc = ra[["Rc"]],
               n_ctrl = sum(ctrl), n_R = sum(reg_c$R), n_Rc = sum(reg_c$Rc),
               n_R_far = sum(reg_c$R_far), n_Rc_far = sum(reg_c$Rc_far))
  })
  do.call(rbind, rows)
}
cat("(A) surface RMSE from stored surfaces, R =", N_ALL, "\n")
surf <- do.call(rbind, parallel::mclapply(seq_len(N_ALL), function(r) {
  do.call(rbind, lapply(SCEN, surface_one, r = r))
}, mc.cores = CORES))
surf <- surf[order(surf$scenario, surf$method, surf$rep), ]
write.csv(surf, file.path(OUT, "surface_per_replicate.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# (B) pointwise ESS from refits of LRC-BART-100
# ---------------------------------------------------------------------------
refit_one <- function(sc, r) {
  out_f <- file.path(OUT, sc, sprintf("rep_%04d.rds", r))
  if (file.exists(out_f)) return(readRDS(out_f))
  t0 <- Sys.time()
  dat <- gen_data_id(sc, "gaussian", seed = rep_seed(SEED, r), n_mc = 0)
  set.seed(fit_seed(SEED, r, METHOD))
  p <- lrc_refit_parts(dat)
  cal <- p$cal; ctrl <- p$ctrl; xc <- p$xc
  tau0 <- tau0_prior_draws(cal)
  pr <- pointwise_prior_ess(cal, xc, tau0)
  re <- pointwise_realized_ess(ctrl, cal, xc, pr$V_f)
  ess_real <- lrc_ess_realized(ctrl, dat$x_rct)
  map <- lrc_borrowing_map(ctrl, xc, type = "abs_g", c = 0.5)$mean
  # check against the stored full-study fit
  stored <- readRDS(result_path(FULL, sc, METHOD, r))
  surf_refit <- colMeans(predict(ctrl, dat$x_rct, arm = "rct_control"))
  x5 <- xc[, 5]
  x5_all <- dat$x_rct[, 5]; x_grid <- rbind(xmat(dat$rwd_ctrl), dat$x_rct)
  reg <- rbind(cbind(region = "R", regional_ess(cal, ctrl, dat$x_rct[x5_all > 2, , drop = FALSE], x_grid, tau0)),
               cbind(region = "Rc", regional_ess(cal, ctrl, dat$x_rct[x5_all <= 2, , drop = FALSE], x_grid, tau0)),
               cbind(region = "all", regional_ess(cal, ctrl, dat$x_rct, x_grid, tau0)))
  reg <- data.frame(scenario = sc, rep = r, reg)
  prof <- data.frame(scenario = sc, rep = r, x5 = x5, inR = x5 > 2, pr, re, map = map)
  res <- list(scenario = sc, rep = r, profiles = prof, regional = reg,
              cal = list(s0_sq = cal$s0_sq, sigma1_sq = cal$sigma1_sq, c_g = cal$c_g,
                         V_mu_f = cal$V_mu_f, ess0_at_s0 = cal$ess0_at_s0,
                         ess0_whole = ess0_whole_chain(cal, tau0), ceiling = cal$ceiling,
                         ess_capped = cal$ess_capped, N_target = cal$N_target,
                         ess_realized = ess_real$mean, ess_realized_cg = ess_real$mean_cg_approx,
                         tau0_sq_post = mean(ctrl$tau0_sq), tau1_sq = mean(ctrl$tau1_sq),
                         w_post = mean(ctrl$w), sigma1_sq_post = mean(ctrl$sigma1_sq)),
              check = list(surface_maxabs = max(abs(surf_refit - stored$surface)),
                           map_maxabs = max(abs(map - stored$map)),
                           ess_prior_stored = stored$ess$prior, ess_realized_stored = stored$ess$realized),
              seconds = as.numeric(Sys.time() - t0, units = "secs"))
  dir.create(dirname(out_f), showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, out_f)
  res
}
cat("(B) refits, R =", N_REP, "on", CORES, "cores\n")
t0 <- Sys.time()
jobs <- expand.grid(r = seq_len(N_REP), sc = SCEN, stringsAsFactors = FALSE)
fits <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  tryCatch(refit_one(jobs$sc[j], jobs$r[j]), error = function(e) {
    cat("ERROR", jobs$sc[j], jobs$r[j], conditionMessage(e), "\n"); NULL })
}, mc.cores = CORES, mc.preschedule = FALSE)
wall_refit <- as.numeric(Sys.time() - t0, units = "mins")
cat("refits done in", round(wall_refit, 1), "min;", sum(!vapply(fits, is.null, TRUE)), "of", nrow(jobs), "\n")
fits <- Filter(Negate(is.null), fits)

prof <- do.call(rbind, lapply(fits, `[[`, "profiles"))
write.csv(prof, file.path(OUT, "pointwise_per_profile.csv"), row.names = FALSE)
calib <- do.call(rbind, lapply(fits, function(f) data.frame(scenario = f$scenario, rep = f$rep, as.data.frame(f$cal),
                                                            as.data.frame(f$check), seconds = f$seconds)))
write.csv(calib, file.path(OUT, "calibration_per_replicate.csv"), row.names = FALSE)
regional <- do.call(rbind, lapply(fits, `[[`, "regional"))
write.csv(regional, file.path(OUT, "regional_ess_per_replicate.csv"), row.names = FALSE)
sumR <- aggregate(cbind(n_S, V_S, c_g_S, ess0_S, ceiling_S, ess_real_S) ~ scenario + region, regional, mean)
sumR$ess0_S_median <- aggregate(ess0_S ~ scenario + region, regional, median)$ess0_S
sumR$ess_real_S_median <- aggregate(ess_real_S ~ scenario + region, regional, median)$ess_real_S
write.csv(sumR, file.path(HERE, "regional_ess_summary.csv"), row.names = FALSE)

# per-replicate region summaries of the pointwise ESS
ess_rep <- do.call(rbind, lapply(fits, function(f) {
  p <- f$profiles
  one <- function(sel, lab) data.frame(scenario = f$scenario, rep = f$rep, region = lab, n = sum(sel),
                                       ess0_mean = mean(p$ess0_x[sel]), ess0_median = median(p$ess0_x[sel]),
                                       ceiling_mean = mean(p$ceiling_x[sel]),
                                       ess_real_mean = mean(p$ess_real_x[sel]), ess_real_median = median(p$ess_real_x[sel]),
                                       ess_post_mean = mean(p$ess_post_x[sel]), ess_post_median = median(p$ess_post_x[sel]),
                                       map_mean = mean(p$map[sel]), g_mean = mean(p$g_mean_x[sel]))
  rbind(one(p$inR, "R"), one(!p$inR, "Rc"), one(p$x5 > 2 + FAR, "R_far"), one(p$x5 < 2 - FAR, "Rc_far"), one(rep(TRUE, nrow(p)), "all"))
}))
write.csv(ess_rep, file.path(OUT, "pointwise_ess_per_replicate.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------
se <- function(v) sd(v) / sqrt(length(v))
summarize_surface <- function(d) {
  agg <- aggregate(cbind(rmse_all, rmse_R, rmse_Rc, rmse_R_far, rmse_Rc_far,
                         rmse_allprof_all, rmse_allprof_R, rmse_allprof_Rc,
                         n_ctrl, n_R, n_Rc, n_R_far, n_Rc_far) ~ scenario + method, d, mean)
  agg$n_rep <- aggregate(rep ~ scenario + method, d, length)$rep
  agg
}
paired_surface <- function(d) {
  out <- list()
  for (sc in unique(d$scenario)) for (reg in c("rmse_all", "rmse_R", "rmse_Rc", "rmse_R_far", "rmse_Rc_far", "rmse_allprof_R", "rmse_allprof_Rc")) {
    a <- d[d$scenario == sc & d$method == "LRC-BART-100", ]; b <- d[d$scenario == sc & d$method == "BART-PP", ]
    b <- b[match(a$rep, b$rep), ]
    dd <- a[[reg]] - b[[reg]]
    out[[length(out) + 1]] <- data.frame(scenario = sc, metric = reg, n_rep = length(dd),
                                         lrc = mean(a[[reg]]), pp = mean(b[[reg]]),
                                         diff = mean(dd), se_paired = se(dd), z = mean(dd) / se(dd),
                                         frac_lrc_better = mean(dd < 0))
  }
  do.call(rbind, out)
}
surf200 <- surf[surf$rep <= N_REP, ]
sumA <- summarize_surface(surf200)
sumA500 <- summarize_surface(surf)
paired <- paired_surface(surf200)
paired500 <- paired_surface(surf)

# pointwise ESS summary (LRC-BART-100): mean over replicates of the per-replicate region means/medians
sumB <- aggregate(cbind(n, ess0_mean, ess0_median, ceiling_mean, ess_real_mean, ess_real_median,
                        ess_post_mean, ess_post_median, map_mean, g_mean) ~ scenario + region, ess_rep, mean)
sumB$n_rep <- aggregate(rep ~ scenario + region, ess_rep, length)$rep
sumC <- aggregate(cbind(s0_sq, sigma1_sq, c_g, V_mu_f, ess0_at_s0, ess0_whole, ceiling, ess_capped, N_target,
                        ess_realized, tau0_sq_post, w_post, surface_maxabs, map_maxabs, seconds) ~ scenario, calib, mean)
sumC$s0_sq_median <- aggregate(s0_sq ~ scenario, calib, median)$s0_sq
sumC$ess0_at_s0_median <- aggregate(ess0_at_s0 ~ scenario, calib, median)$ess0_at_s0
sumC$surface_maxabs_max <- aggregate(surface_maxabs ~ scenario, calib, max)$surface_maxabs

# summary.csv: scenario x method x region (long)
wide_to_long <- function(a) {
  do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    data.frame(scenario = a$scenario[i], method = a$method[i], n_rep = a$n_rep[i],
               region = c("all", "R", "Rc", "R_far", "Rc_far"),
               n_profiles = unlist(a[i, c("n_ctrl", "n_R", "n_Rc", "n_R_far", "n_Rc_far")]),
               surf_rmse = unlist(a[i, c("rmse_all", "rmse_R", "rmse_Rc", "rmse_R_far", "rmse_Rc_far")]))
  }))
}
summary_long <- wide_to_long(sumA)
key <- paste(summary_long$scenario, summary_long$region)
keyB <- paste(sumB$scenario, sumB$region)
for (v in c("ess0_mean", "ess0_median", "ceiling_mean", "ess_real_mean", "ess_real_median", "ess_post_mean", "ess_post_median", "map_mean", "g_mean")) {
  summary_long[[paste0("lrc_", v)]] <- ifelse(summary_long$method == "LRC-BART-100", sumB[[v]][match(key, keyB)], NA)
}
rownames(summary_long) <- NULL
write.csv(summary_long, file.path(HERE, "summary.csv"), row.names = FALSE)
write.csv(wide_to_long(sumA500), file.path(HERE, "summary_R500.csv"), row.names = FALSE)
write.csv(paired, file.path(HERE, "paired.csv"), row.names = FALSE)
write.csv(paired500, file.path(HERE, "paired_R500.csv"), row.names = FALSE)
write.csv(sumC, file.path(HERE, "calibration_summary.csv"), row.names = FALSE)
wall_total <- as.numeric(Sys.time() - t_start, units = "mins")
writeLines(c(sprintf("start %s", format(t_start)), sprintf("N_REP %d CORES %d", N_REP, CORES),
             sprintf("wall_minutes_refits %.1f", wall_refit), sprintf("wall_minutes_total %.1f", wall_total),
             sprintf("mean_seconds_per_refit %.2f", mean(calib$seconds)),
             sprintf("max_surface_maxabs %.3g", max(calib$surface_maxabs)),
             sprintf("max_map_maxabs %.3g", max(calib$map_maxabs))),
           file.path(HERE, "wall_time.txt"))
options(width = 200)
cat("\n== surface RMSE at RCT control profiles (R =", N_REP, ") ==\n"); print(sumA, digits = 3)
cat("\n== paired LRC-BART-100 minus BART-PP ==\n"); print(paired, digits = 3)
cat("\n== pointwise ESS (LRC-BART-100) ==\n"); print(sumB, digits = 3)
cat("\n== regional ESS (LRC-BART-100) ==\n"); print(sumR, digits = 3)
cat("\n== calibration ==\n"); print(sumC, digits = 3)
cat("\nwall", round(wall_total, 1), "min\n")
