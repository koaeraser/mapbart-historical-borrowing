# compare_cahb.R: pair the CAHB rows (results/per_replicate.csv, R = 200,
# seed 2026) with replicates 1..200 of the full Gaussian study
# (../full/gaussian/per_replicate.csv, R = 500, seed 2026; identical data
# seeds 2026 + r) and write
#   summary_table.csv  one row per scenario x method: ATE bias, SD, RMSE (MC
#                      RMSE of the posterior mean), coverage, power (95% CrI
#                      excludes 0), control-surface RMSE overall / in R / in
#                      R^c, map summaries, ESS, seconds; n_rep = 200
#   paired.csv         per scenario x method, differences from LRC-BART-100
#                      over the same 200 replicates with paired MC standard
#                      errors: MSE difference (mean of the per-replicate
#                      squared-error difference, SE = sd / sqrt(R)), the RMSE
#                      difference with its delta-method SE, and the bias, SD,
#                      coverage, power and surface-RMSE differences.
# The pairing is verified by requiring identical ate_truth per (scenario,
# rep) in both files.
# Usage: Rscript compare_cahb.R
suppressPackageStartupMessages(library(data.table))
HERE <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
FULL <- file.path(HERE, "..", "full", "gaussian", "per_replicate.csv")
scen <- c("Sc1", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d1", "Sc5_d2")
ref_methods <- c("LRC-BART-100", "BART-PP", "BART-CP", "BART-NP", "PSCL", "MAP-100")
cols <- c("scenario", "method", "rep", "ate_pm", "ate_truth", "ate_sd", "ate_cover", "ate_reject", "ate_reject90",
          "ctrl_pm", "ctrl_truth", "surf_rmse", "surf_rmse_R", "surf_rmse_Rc",
          "map_R", "map_Rc", "map_all", "ess_prior", "ess_realized", "seconds")
full <- fread(FULL, select = cols)[scenario %in% scen & method %in% ref_methods & rep <= 200]
cahb <- fread(file.path(HERE, "results", "per_replicate.csv"), select = cols)[scenario %in% scen & rep <= 200]
stopifnot(nrow(full) == length(scen) * length(ref_methods) * 200)
# pairing check: same truth per (scenario, rep)
tr_full <- unique(full[, .(scenario, rep, ate_truth)])
tr_cahb <- unique(cahb[, .(scenario, rep, ate_truth)])
chk <- merge(tr_full, tr_cahb, by = c("scenario", "rep"), suffixes = c("_full", "_cahb"))
stopifnot(nrow(chk) == length(scen) * 200, max(abs(chk$ate_truth_full - chk$ate_truth_cahb)) < 1e-8)
cat("pairing check passed: identical ate_truth on", nrow(chk), "(scenario, rep) pairs\n")
all <- rbind(full, cahb)
all[, err := ate_pm - ate_truth]
all[, sq := err^2]

summ <- all[, .(n_rep = .N, bias = mean(err), sd = mean(ate_sd), rmse = sqrt(mean(sq)),
                cov = mean(ate_cover), pow = mean(ate_reject), pow90 = mean(ate_reject90),
                ctrl_bias = mean(ctrl_pm - ctrl_truth), ctrl_rmse = sqrt(mean((ctrl_pm - ctrl_truth)^2)),
                surf_rmse = mean(surf_rmse), surf_rmse_R = mean(surf_rmse_R, na.rm = TRUE),
                surf_rmse_Rc = mean(surf_rmse_Rc, na.rm = TRUE),
                map_R = mean(map_R, na.rm = TRUE), map_Rc = mean(map_Rc, na.rm = TRUE), map_all = mean(map_all, na.rm = TRUE),
                ess_prior = mean(ess_prior, na.rm = TRUE), ess_realized = mean(ess_realized, na.rm = TRUE),
                seconds = mean(seconds)),
            by = .(scenario, method)]
lev <- c("CAHB", "CAHB-p10", "CAHB-lit", ref_methods)
summ[, method := factor(method, levels = lev)]; setorder(summ, scenario, method); summ[, method := as.character(method)]
fwrite(summ, file.path(HERE, "summary_table.csv"))

# paired differences from LRC-BART-100
L <- all[method == "LRC-BART-100", .(scenario, rep, sq_L = sq, err_L = err, sd_L = ate_sd, cov_L = ate_cover,
                                      rej_L = ate_reject, surf_L = surf_rmse, surfR_L = surf_rmse_R, surfRc_L = surf_rmse_Rc)]
P <- merge(all[method != "LRC-BART-100"], L, by = c("scenario", "rep"))
se <- function(x) sd(x) / sqrt(length(x))
paired <- P[, {
  mse_m <- mean(sq); mse_L <- mean(sq_L); d <- sq - sq_L
  rmse_m <- sqrt(mse_m); rmse_L <- sqrt(mse_L)
  list(n_rep = .N, rmse = rmse_m, rmse_lrc = rmse_L,
       mse_diff = mse_m - mse_L, mse_diff_se = se(d),
       rmse_diff = rmse_m - rmse_L, rmse_diff_se = se(d) / (rmse_m + rmse_L),
       bias_diff = mean(err - err_L), bias_diff_se = se(err - err_L),
       sd_diff = mean(ate_sd - sd_L), sd_diff_se = se(ate_sd - sd_L),
       cov_diff = mean(ate_cover - cov_L), cov_diff_se = se(ate_cover - cov_L),
       pow_diff = mean(ate_reject - rej_L), pow_diff_se = se(ate_reject - rej_L),
       surf_diff = mean(surf_rmse - surf_L), surf_diff_se = se(surf_rmse - surf_L),
       surfR_diff = mean(surf_rmse_R - surfR_L, na.rm = TRUE), surfR_diff_se = se(na.omit(surf_rmse_R - surfR_L)),
       surfRc_diff = mean(surf_rmse_Rc - surfRc_L, na.rm = TRUE), surfRc_diff_se = se(na.omit(surf_rmse_Rc - surfRc_L)))
}, by = .(scenario, method)]
paired[, method := factor(method, levels = lev)]; setorder(paired, scenario, method); paired[, method := as.character(method)]
fwrite(paired, file.path(HERE, "paired.csv"))
cat("wrote summary_table.csv (", nrow(summ), "rows ) and paired.csv (", nrow(paired), "rows )\n")
print(summ[, .(scenario, method, bias = round(bias, 3), sd = round(sd, 3), rmse = round(rmse, 3), cov = round(cov, 3),
               pow = round(pow, 3), surf = round(surf_rmse, 3), surf_R = round(surf_rmse_R, 3), surf_Rc = round(surf_rmse_Rc, 3),
               map_R = round(map_R, 3), map_Rc = round(map_Rc, 3), ess_r = round(ess_realized, 1))], nrows = 100)
