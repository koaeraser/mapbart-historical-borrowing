source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
suppressPackageStartupMessages({library(BART); library(survival)})
CODE_ROOT <- "/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/MAP-BART/revision/01-code"
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
R <- 20
# ---- repo post-processing (BARTv1.R survival lines 325-361, 470-560), verbatim logic ----
S_mix_fun <- function(t, mu_vec, sig, w) sum(w * pnorm(log(t), mean = mu_vec, sd = sig, lower.tail = FALSE))
repo_pop_median <- function(mu_draws, sig_draw, lower) {
  w <- rep(1/ncol(mu_draws), ncol(mu_draws)); out <- rep(NA_real_, nrow(mu_draws))
  upper_default <- pmin(exp(apply(mu_draws, 1, max) + 6 * sig_draw), 1e10)
  for (m in seq_len(nrow(mu_draws))) {
    f <- function(t) S_mix_fun(t, mu_draws[m, ], sig_draw[m], w) - 0.5
    if (f(upper_default[m]) > 0) { out[m] <- NA; next }
    lo <- max(1e-10, lower)
    out[m] <- if (f(lo) < 0) lo else uniroot(f, interval = c(lo, upper_default[m]))$root
  }
  out
}
# repo BART settings (BARTv1.R): ntree 50, ndpost 1500, nskip 100, k 2, power 2, base 0.95, sigdf 3, sigquant 0.9,
# sigest from a lognormal survreg on the training data; the C++ (cabart.cpp) is not in the repo, BART::abart stands in.
repo_fit <- function(X, times, delta, X_test, seed) {
  aft <- survreg(Surv(times, delta) ~ ., data = data.frame(X), dist = "lognormal")
  set.seed(seed)
  fit <- abart(x.train = X, times = times, delta = delta, x.test = X_test, K = 2, ntree = 50L, k = 2, power = 2, base = 0.95,
               sigdf = 3, sigquant = 0.9, sigest = aft$scale, ndpost = 1500L, nskip = 100L, keepevery = 1L, printevery = 100000L)
  list(mu = fit$yhat.test, sigma = tail(fit$sigma, 1500L))
}
summ <- function(ps, eff) { ps <- ps[is.finite(ps)]
  c(med = median(ps) - eff, mean = mean(ps) - eff, sd = sd(ps), lo = quantile(ps, .025), hi = quantile(ps, .975),
    cov = as.numeric(quantile(ps, .025) <= eff & quantile(ps, .975) >= eff), pow = as.numeric(mean(ps > 1) > 0.95)) }
res_repo <- res_ours <- NULL
for (r in 1:R) {
  # ---- repo DGP replicate r, repo pipeline ----
  d <- readRDS(file.path(SCR, "mapbart-sim-survival/data_v2", sprintf("data_p10_sc1E_alternative_%d.RData", r)))
  X <- as.matrix(d$X[, paste0("X", 1:10)]); D <- d$X$D; Z <- d$X$Z
  ic <- Z == 0 & D == 1; it <- Z == 1; itest <- D == 1
  fc <- repo_fit(X[ic, ], d$y[ic], d$delta[ic], X[itest, ], 1000 + r)
  ft <- repo_fit(X[it, ], d$y[it], d$delta[it], X[itest, ], 2000 + r)
  lower <- min(d$y[d$delta == 1]) / 10
  ps <- repo_pop_median(ft$mu, ft$sigma, lower) / repo_pop_median(fc$mu, fc$sigma, lower)
  eff <- exp(d$treat_eff_true)
  # our estimand code applied to the same draws
  ps_ours_code <- standardized_draws(ft$mu, fc$mu, "survival", ft$sigma, fc$sigma, 3)$ate_med
  res_repo <- rbind(res_repo, c(rep = r, n_ctrl = sum(ic), eff = eff, summ(ps, eff), ours_code_mean = mean(ps_ours_code) - eff,
                                ours_code_med = median(ps_ours_code) - eff))
  # ---- our DGP replicate r (seed 2026 + r, as in the full study), our harness BART-NP ----
  dat <- gen_data(1, "survival", seed = rep_seed(2026, r), n_mc = 0)
  set.seed(fit_seed(2026, r, "BART-NP")); fit <- fit_bart_np(dat, "survival")
  ps2 <- fit$draws$ate_med; eff2 <- dat$truth$ate_med
  res_ours <- rbind(res_ours, c(rep = r, n_ctrl = nrow(dat$rct_ctrl), eff = eff2, summ(ps2, eff2)))
  cat(sprintf("rep %2d repo: med-bias %.3f mean-bias %.3f sd %.2f | ours: med-bias %.3f mean-bias %.3f sd %.2f\n", r,
              res_repo[r, "med"], res_repo[r, "mean"], res_repo[r, "sd"], res_ours[r, "med"], res_ours[r, "mean"], res_ours[r, "sd"]))
}
agg <- function(m) c(bias_median = mean(m[, "med"]), bias_mean = mean(m[, "mean"]), sd = mean(m[, "sd"]),
                     rmse_median = sqrt(mean(m[, "med"]^2)), rmse_mean = sqrt(mean(m[, "mean"]^2)), cov = mean(m[, "cov"]), pow = mean(m[, "pow"]))
out <- rbind(repo_dgp_repo_pipeline = agg(res_repo), ours_dgp_our_harness = agg(res_ours))
print(round(out, 3))
cat("repo draws through our estimands.R: mean-bias", round(mean(res_repo[, "ours_code_mean"]), 3), " median-bias", round(mean(res_repo[, "ours_code_med"]), 3), "\n")
cat("mean n_ctrl repo", mean(res_repo[, "n_ctrl"]), " ours", mean(res_ours[, "n_ctrl"]), "\n")
saveRDS(list(repo = res_repo, ours = res_ours, agg = out), file.path(SCR, "out/q3_bart_np.rds"))
