# replay_extras.R: re-fit CAHB (primary) on the study replicates with the
# harness seeds to recover the per-replicate diagnostics that fit_one() does
# not store (kernel covariates, bandwidth multiplier, phi_hat^2, lambda_1,
# lambda_2, gamma^2, iterations). Writes extras_cahb.csv.
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
HERE <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT); source(file.path(CODE_ROOT, "R", "cahb.R"))
scen <- c("Sc1", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d1", "Sc5_d2")
jobs <- expand.grid(r = 1:200, scenario = scen, stringsAsFactors = FALSE)
one <- function(j) {
  sc <- jobs$scenario[j]; r <- jobs$r[j]
  dat <- gen_data_id(sc, "gaussian", seed = rep_seed(2026, r), n_mc = 0)
  set.seed(fit_seed(2026, r, "CAHB"))
  f <- fit_cahb(dat, "gaussian"); e <- f$extra
  data.frame(scenario = sc, rep = r, vars = paste0("X", e$vars, collapse = "+"), mult = e$mult_c,
             sw_med = median(e$sum_w_ctrl), sw_min = min(e$sum_w_ctrl), phi2_hat = e$phi2_hat, phi0 = e$phi0,
             gamma2 = e$gamma2, lam1 = e$lam1, lam2 = e$lam2, iter = e$iter, n_tau_pos = sum(e$tau > 0),
             ess_realized = f$ess$realized, ate_pm = mean(f$draws$ate))
}
out <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = 12))
data.table::fwrite(out, file.path(HERE, "extras_cahb.csv"))
# consistency with the stored run
pr <- data.table::fread(file.path(HERE, "results", "per_replicate.csv"))[method == "CAHB", .(scenario, rep, ess_stored = ess_realized, ate_stored = ate_pm)]
m <- merge(data.table::as.data.table(out), pr, by = c("scenario", "rep"))
cat("max |ess replay - stored| =", max(abs(m$ess_realized - m$ess_stored)), " max |ate replay - stored| =", max(abs(m$ate_pm - m$ate_stored)), "\n")
suppressPackageStartupMessages(library(data.table)); o <- as.data.table(out)
cat("\nkernel covariate sets (all scenarios):\n"); print(sort(table(o$vars), decreasing = TRUE)[1:8])
cat("\nfrequency each covariate is in the kernel:\n"); print(round(sapply(paste0("X", 1:10), function(v) mean(grepl(paste0("(^|\\+)", v, "($|\\+)"), o$vars))), 3))
print(o[, .(mult = round(mean(mult), 2), sw_med = round(mean(sw_med), 1), sw_min = round(median(sw_min), 2), phi2_hat = round(mean(phi2_hat), 2),
            gamma2 = round(mean(gamma2), 3), lam1 = round(mean(lam1), 2), lam2 = round(mean(lam2), 1), iter = round(mean(iter), 1),
            iter_max = max(iter), frac_tau_pos = round(mean(n_tau_pos / 300), 2)), by = scenario])
