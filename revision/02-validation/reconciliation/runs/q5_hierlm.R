source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
t0 <- Sys.time()
for (sc in c(1, 2)) {
  run_repo_script("mapbart-sim-gaussian-single-arm", "HierLM.R",
    list(seed_method = 789, RCT_ctrl = FALSE, prior_vals = c(0.05, 0.5), sc = sc, subsc = "E", hypo = "alternative"))
  for (pr in c(0.05, 0.5)) {
    h <- readRDS(file.path(SCR, sprintf("mapbart-sim-gaussian-single-arm/res/HierLM_p10_sc%dE_prior%s_alternative.RData", sc, pr)))
    h <- h[!is.na(h$bias), ]
    cat(sprintf("[repo HierLM single-arm Sc%d prior s2=%s -> IG(1.5, %.3f), R=%d] bias %.3f sd %.3f rmse %.3f cov %.2f pow %.2f ctrl bias %.3f ctrl rmse(draw) %.3f\n",
        sc, pr, 1.5 * pr, nrow(h), mean(h$bias), mean(h$sd), sqrt(mean(h$rmse)), mean(h$coverage), mean(h$tp),
        mean(h$bias.ctrl.median.pop, na.rm = TRUE), mean(h$w2distance.ctrl.median.pop, na.rm = TRUE)))
  }
}
cat("elapsed", format(Sys.time() - t0), "\n")
