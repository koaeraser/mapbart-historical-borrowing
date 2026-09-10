source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
t0 <- Sys.time()
# Two-arm survival Sc1, run_all.R overrides (seed_rct 123, seed_rwd 456, fresh RWD), 20 reps.
# Stop before the R^2 testing block at the end of the script.
run_repo_script("mapbart-sim-survival", "data_gen_p10_v2.R",
  list(seed_rct = 123, seed_rwd = 456, rwd_frozen = FALSE, p_obs = 10, niter = 20, hypo = "alternative", sc = 1, subsc = "E"),
  stop_at = "## Testing code")
# Single-arm Gaussian Sc2 and Sc1 (data_gen_p10_v3.R), run_all.R overrides, 20 reps.
for (s in c(2, 1)) run_repo_script("mapbart-sim-gaussian-single-arm", "data_gen_p10_v3.R",
  list(seed_rct = 123, seed_rwd = 456, rwd_frozen = FALSE, p_obs = 10, niter = 20, hypo = "alternative", sc = s, subsc = "E"),
  stop_at = "## Testing code")
cat("files:", length(list.files(file.path(SCR, "mapbart-sim-survival/data_v2"))),
    length(list.files(file.path(SCR, "mapbart-sim-gaussian-single-arm/data_v3"))), "\n")
d <- readRDS(file.path(SCR, "mapbart-sim-survival/data_v2/data_p10_sc1E_alternative_1.RData"))
cat("surv sc1 rep1: n_rct", d$n_rct, "n_rwd", d$n_rwd, "trt", sum(d$X$D==1 & d$X$Z==1), "ctrl", sum(d$X$D==1 & d$X$Z==0),
    "true ratio", exp(d$treat_eff_true), "event rate", mean(d$delta), "\n")
g <- readRDS(file.path(SCR, "mapbart-sim-gaussian-single-arm/data_v3/data_p10_sc2E_alternative_1.RData"))
cat("gauss sa sc2 rep1: n_rct", g$n_rct, "n_rwd", g$n_rwd, "Z all 1:", all(g$X$Z[g$X$D==1]==1), "true ATE", g$treat_eff_true,
    "mean ctrl truth", g$true_mean_ctrl_pop, "mean y rwd", mean(g$y[g$X$D==0]), "\n")
cat("elapsed", format(Sys.time() - t0), "\n")
