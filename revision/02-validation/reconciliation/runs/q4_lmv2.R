source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
suppressPackageStartupMessages(library(BART))
CODE_ROOT <- "/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/MAP-BART/revision/01-code"
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R")); source(file.path(CODE_ROOT, "R", "single_arm.R")); register_single_arm_methods()
R <- 20
# 1. repo LMv2.R verbatim (single-arm, rwd_w = 1 = LM-CP), Sc2, run_all.R seeds
t0 <- Sys.time()
run_repo_script("mapbart-sim-gaussian-single-arm", "LMv2.R",
  list(seed_method = 789, RCT_ctrl = FALSE, rwd_w_vals = 1, sc = 2, subsc = "E", hypo = "alternative"))
lm <- readRDS(file.path(SCR, "mapbart-sim-gaussian-single-arm/res/LMv2_p10_sc2E_alternative_w1.RData"))
lm <- lm[!is.na(lm$bias), ]
cat("\n[repo LMv2 Sc2 single-arm, R =", nrow(lm), "] bias", round(mean(lm$bias), 3), "sd", round(mean(lm$sd), 3),
    "rmse", round(sqrt(mean(lm$rmse)), 3), "cov", round(mean(lm$coverage), 2), "pow", round(mean(lm$tp), 2),
    "ctrl bias", round(mean(lm$bias.ctrl.median.pop, na.rm = TRUE), 3), " elapsed", format(Sys.time() - t0), "\n")
# 2. our LM-CP and BART-CP on the repo's v3 data (converted to the harness data object)
to_dat <- function(d) {
  X <- as.matrix(d$X[, paste0("X", 1:10)]); D <- d$X$D
  df <- as.data.frame(X); df$d <- D; df$z <- d$X$Z; df$y <- d$y
  list(rct_trt = df[D == 1, ], rct_ctrl = df[D == 1 & d$X$Z == 0, ], rwd_ctrl = df[D == 0, ],
       x_rct = X[D == 1, ], z_rct = d$X$Z[D == 1],
       truth = list(ate = d$treat_eff_true, trt = d$true_mean_trt_pop, ctrl = d$true_mean_ctrl_pop),
       mu0_rct = d$true_mean_ctrl, region_rct = rep(FALSE, sum(D == 1)),
       constants = list(outcome = "gaussian", tau = 3, scenario_id = "Sc2"))
}
res <- list()
for (sc in c(2, 1)) for (r in 1:R) {
  d <- readRDS(file.path(SCR, "mapbart-sim-gaussian-single-arm/data_v3", sprintf("data_p10_sc%dE_alternative_%d.RData", sc, r)))
  dat <- to_dat(d)
  for (m in c("LM-CP", "BART-CP")) { set.seed(fit_seed(2026, r, m)); f <- get_method(m)(dat, "gaussian")
    mt <- compute_metrics(f, dat); res[[length(res) + 1]] <- data.frame(sc = sc, method = m, src = "repo_v3_data", rep = r, mt) }
  if (sc == 2) { dat2 <- gen_data(2, "gaussian", n_rct = 200, n_rwd = 300, ratio = .Machine$double.xmax, seed = rep_seed(2026, r), n_mc = 0)
    for (m in c("LM-CP")) { set.seed(fit_seed(2026, r, m)); f <- get_method(m)(dat2, "gaussian")
      mt <- compute_metrics(f, dat2); res[[length(res) + 1]] <- data.frame(sc = sc, method = m, src = "our_dgp", rep = r, mt) } }
}
res <- do.call(rbind, res)
ag <- aggregate(cbind(ate_bias, ate_sd, ate_cover, ctrl_bias) ~ sc + method + src, res, mean)
ag$ate_rmse <- aggregate(ate_bias ~ sc + method + src, res, function(b) sqrt(mean(b^2)))$ate_bias
print(ag, digits = 3)
saveRDS(list(repo_lmv2 = lm, ours = res, agg = ag), file.path(SCR, "out/q4_lmv2.rds"))
