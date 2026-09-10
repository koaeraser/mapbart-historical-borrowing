# Phase 2b tuning of the discrepancy ensemble g (STATE.md, 2026-09-09).
# Runs named LRC-BART configurations on one or more scenarios through the
# harness (paired replicate data, BART-PP as the benchmark in every run) and
# writes a per-config summary with paired differences against BART-PP.
#
# Usage: Rscript run_factorial.R <stage> [R] [workers]
#   stage "ofat"    one-factor-at-a-time from the baseline on Sc4_X5_d2
#   stage "combo"   the promising combinations on Sc4_X5_d2
#   stage "confirm" the chosen configuration on Sc4_X5_d1, Sc3_rho0, Sc1, Sc5_d2
# Config names encode the factors: H<H_g>_a<alpha_g>b<beta_g>_t<tau1_mult>_n<n_min_g>
# with optional suffixes _warm (warm start from the BART-PP offset) and
# _sw<k> (k g-sweeps per iteration).
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
HERE <- file.path(dirname(CODE_ROOT), "02-validation", "tuning_phase2b")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "ofat"
R <- if (length(args) > 1) as.integer(args[2]) else 40L
workers <- if (length(args) > 2) as.integer(args[3]) else 16L

cfg <- function(H_g = 10, alpha_g = 0.5, beta_g = 3, tau1_mult = 1, n_min_g = 5,
                warm_start = FALSE, g_sweeps = 1, N_target = 100, w_prior = c(4, 1)) {
  nm <- sprintf("H%d_a%gb%g_t%g_n%d", H_g, alpha_g, beta_g, tau1_mult, n_min_g)
  if (warm_start) nm <- paste0(nm, "_warm")
  if (g_sweeps > 1) nm <- paste0(nm, "_sw", g_sweeps)
  if (!identical(w_prior, c(4, 1))) nm <- paste0(nm, "_w", w_prior[1], w_prior[2])
  if (!identical(N_target, 100)) nm <- paste0(nm, "_N", paste(N_target, collapse = ""))
  list(name = nm, args = list(H_g = H_g, alpha_g = alpha_g, beta_g = beta_g, tau1_mult = tau1_mult,
                              n_min_g = n_min_g, warm_start = warm_start, g_sweeps = g_sweeps,
                              N_target = N_target, w_prior = w_prior))
}
register_cfg <- function(c) {
  a <- c$args
  register_method(c$name, function(dat, outcome, ...) do.call(fit_lrcbart, c(list(dat = dat, outcome = outcome), a)))
}

configs <- switch(stage,
  ofat = list(cfg(),                                   # baseline (dev-run configuration)
              cfg(H_g = 20), cfg(H_g = 30),
              cfg(alpha_g = 0.95, beta_g = 2), cfg(alpha_g = 0.95, beta_g = 1),
              cfg(tau1_mult = 4), cfg(tau1_mult = 16),
              cfg(n_min_g = 1),
              cfg(warm_start = TRUE),                  # mixing diagnostic
              cfg(g_sweeps = 5)),                      # mixing fix candidate
  # After the OFAT stage (none of H_g in {10, 20, 30}, the split prior, the
  # slab scale, n_min, a warm start or extra sweeps moved g) and the oracle
  # runs (g under-shoots even with f known, and more trees fuzz the map), the
  # combinations tested are fewer g trees with the other factors at their
  # promising levels.
  combo = list(cfg(H_g = 1), cfg(H_g = 2), cfg(H_g = 3), cfg(H_g = 5),
               cfg(H_g = 3, alpha_g = 0.95, beta_g = 2), cfg(H_g = 5, alpha_g = 0.95, beta_g = 2),
               cfg(H_g = 2, alpha_g = 0.95, beta_g = 2),
               cfg(H_g = 3, n_min_g = 1), cfg(H_g = 3, tau1_mult = 4),
               cfg(H_g = 5, tau1_mult = 4), cfg(H_g = 1, alpha_g = 0.95, beta_g = 2)),
  confirm = {
    src <- file.path(HERE, "chosen_config.R"); stopifnot(file.exists(src)); source(src)
    c(list(cfg()), lapply(CANDIDATES, function(a) do.call(cfg, a)))
  },
  # moderate shifts (Sc3, Sc4 delta 1): the slab's Occam cost log(tau1/tau0)
  # + log(w/(1-w)) against the route where f absorbs the RCT; test a cheaper
  # slab (tau1_mult 0.5) and a flatter w prior at H_g 5 and 10
  moderate = list(cfg(), cfg(H_g = 5),
                  cfg(tau1_mult = 0.5), cfg(H_g = 5, tau1_mult = 0.5),
                  cfg(w_prior = c(1, 1)), cfg(H_g = 5, w_prior = c(1, 1)),
                  cfg(H_g = 5, tau1_mult = 0.5, w_prior = c(1, 1))),
  stop("unknown stage"))
scenarios <- switch(stage, ofat = "Sc4_X5_d2", combo = "Sc4_X5_d2",
                    confirm = c("Sc4_X5_d1", "Sc4_X5_d2", "Sc3_rho0", "Sc1", "Sc5_d2"),
                    moderate = c("Sc4_X5_d1", "Sc3_rho0", "Sc4_X5_d2"))
for (c in configs) register_cfg(c)
methods <- c("BART-PP", vapply(configs, `[[`, "", "name"))
out_dir <- file.path(HERE, paste0("results_", stage))
t0 <- Sys.time()
cat("stage", stage, "R", R, "workers", workers, "start", format(t0), "\n")
run_sim(scenarios, "gaussian", methods, R = R, seed = 2026, workers = workers, out_dir = out_dir)
cat("fits done, wall minutes:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
res <- summarize_sim(out_dir)

# Paired summary against BART-PP (same replicate data), one row per scenario x config.
suppressPackageStartupMessages(library(data.table))
pr <- as.data.table(res$per_replicate)
pp <- pr[method == "BART-PP", .(scenario, rep, pp_rmse_R = surf_rmse_R, pp_rmse_Rc = surf_rmse_Rc,
                                pp_rmse = surf_rmse, pp_err = ate_pm - ate_truth)]
d <- merge(pr[method != "BART-PP"], pp, by = c("scenario", "rep"))
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
summ <- d[, .(n = .N,
              g_R = mean(g_mean_R), g_R_q10 = q(g_mean_R, .1), g_R_q90 = q(g_mean_R, .9),
              g_Rc = mean(g_mean_Rc), g_Rc_q10 = q(g_mean_Rc, .1), g_Rc_q90 = q(g_mean_Rc, .9),
              map_R = mean(map_R), map_Rc = mean(map_Rc), map_all = mean(map_all),
              spike_R = mean(mapspike_R), spike_Rc = mean(mapspike_Rc),
              rmse_R = mean(surf_rmse_R), dR = mean(surf_rmse_R - pp_rmse_R), dR_se = sd(surf_rmse_R - pp_rmse_R) / sqrt(.N),
              rmse_Rc = mean(surf_rmse_Rc), dRc = mean(surf_rmse_Rc - pp_rmse_Rc), dRc_se = sd(surf_rmse_Rc - pp_rmse_Rc) / sqrt(.N),
              rmse_surf = mean(surf_rmse), dS = mean(surf_rmse - pp_rmse), dS_se = sd(surf_rmse - pp_rmse) / sqrt(.N),
              bias = mean(ate_pm - ate_truth), bias_se = sd(ate_pm - ate_truth) / sqrt(.N),
              rmse = sqrt(mean((ate_pm - ate_truth)^2)), rmse_pp = sqrt(mean(pp_err^2)),
              dsq = mean((ate_pm - ate_truth)^2 - pp_err^2), dsq_se = sd((ate_pm - ate_truth)^2 - pp_err^2) / sqrt(.N),
              cover = mean(ate_cover), sd = mean(ate_sd),
              ess0 = mean(ess_prior), ess_real = mean(ess_realized), ceiling = mean(ess_ceiling), capped = mean(ess_capped),
              secs = mean(seconds)),
          by = .(scenario, method)]
fwrite(summ, file.path(HERE, paste0("summary_", stage, ".csv")))
f <- function(x, k = 3) formatC(x, format = "f", digits = k)
cat("\n", stage, ": paired against BART-PP (dR = surface RMSE in R minus BART-PP's, SE in parentheses)\n", sep = "")
for (sc in unique(summ$scenario)) {
  cat("\n## ", sc, "\n", sep = "")
  cat("| config | g_R (q10,q90) | g_Rc (q10,q90) | map_R | map_Rc | rmse_R | dR (SE) | rmse_Rc | dRc (SE) | bias (SE) | rmse | rmse_PP | cov | ESSreal | secs |\n")
  cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
  for (i in which(summ$scenario == sc)) {
    r <- summ[i]
    cat(sprintf("| %s | %s (%s, %s) | %s (%s, %s) | %s | %s | %s | %s (%s) | %s | %s (%s) | %s (%s) | %s | %s | %s | %s | %s |\n",
                r$method, f(r$g_R, 2), f(r$g_R_q10, 2), f(r$g_R_q90, 2), f(r$g_Rc, 2), f(r$g_Rc_q10, 2), f(r$g_Rc_q90, 2),
                f(r$map_R, 2), f(r$map_Rc, 2), f(r$rmse_R), f(r$dR), f(r$dR_se), f(r$rmse_Rc), f(r$dRc), f(r$dRc_se),
                f(r$bias), f(r$bias_se), f(r$rmse), f(r$rmse_pp), f(r$cover, 2), f(r$ess_real, 1), f(r$secs, 1)))
  }
}
cat("\nwall minutes total:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "\n")
