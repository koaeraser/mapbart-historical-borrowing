source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
SCR <- sub("/recon$", "/recon100", SCR)
run_repo_script("mapbart-sim-gaussian-single-arm", "data_gen_p10_v3.R",
  list(seed_rct = 123, seed_rwd = 456, rwd_frozen = FALSE, p_obs = 10, niter = 100, hypo = "alternative", sc = 2, subsc = "E"),
  stop_at = "## Testing code")
run_repo_script("mapbart-sim-gaussian-single-arm", "LMv2.R",
  list(seed_method = 789, RCT_ctrl = FALSE, rwd_w_vals = 1, sc = 2, subsc = "E", hypo = "alternative"))
lm <- readRDS(file.path(SCR, "mapbart-sim-gaussian-single-arm/res/LMv2_p10_sc2E_alternative_w1.RData")); lm <- lm[!is.na(lm$bias), ]
cat("\n[repo LMv2 Sc2 single-arm, R =", nrow(lm), "] bias", round(mean(lm$bias), 3), "(MC se", round(sd(lm$bias)/sqrt(nrow(lm)), 3), ") sd", round(mean(lm$sd), 3),
    "rmse", round(sqrt(mean(lm$rmse)), 3), "cov", round(mean(lm$coverage), 2), "pow", round(mean(lm$tp), 2),
    "ctrl bias", round(mean(lm$bias.ctrl.median.pop, na.rm = TRUE), 3), "ctrl rmse(draw)", round(mean(lm$w2distance.ctrl.median.pop, na.rm = TRUE), 3), "\n")
