source("/private/tmp/claude-501/-Users-yuanj/692eaa23-bac6-4a08-87f5-574952eb8c31/scratchpad/recon/helpers.R")
SCR <- sub("/recon$", "/recon100", SCR)
suppressPackageStartupMessages({library(BART); library(survival); library(parallel)})
CODE_ROOT <- "/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/MAP-BART/revision/01-code"
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
R <- 50
run_repo_script("mapbart-sim-survival", "data_gen_p10_v2.R",
  list(seed_rct = 123, seed_rwd = 456, rwd_frozen = FALSE, p_obs = 10, niter = R, hypo = "alternative", sc = 1, subsc = "E"),
  stop_at = "## Testing code")
S_mix_fun <- function(t, mu_vec, sig, w) sum(w * pnorm(log(t), mean = mu_vec, sd = sig, lower.tail = FALSE))
repo_pop_median <- function(mu_draws, sig_draw, lower) {
  w <- rep(1/ncol(mu_draws), ncol(mu_draws)); out <- rep(NA_real_, nrow(mu_draws))
  upper_default <- pmin(exp(apply(mu_draws, 1, max) + 6 * sig_draw), 1e10)
  for (m in seq_len(nrow(mu_draws))) { f <- function(t) S_mix_fun(t, mu_draws[m, ], sig_draw[m], w) - 0.5
    if (f(upper_default[m]) > 0) { out[m] <- NA; next }; lo <- max(1e-10, lower)
    out[m] <- if (f(lo) < 0) lo else uniroot(f, interval = c(lo, upper_default[m]))$root }
  out }
fit_arm <- function(X, times, delta, X_test, seed, sigest_mode, ndpost = 1500L, nskip = 100L) {
  sigest <- if (sigest_mode == "survreg") survreg(Surv(times, delta) ~ ., data = data.frame(X), dist = "lognormal")$scale else NA
  set.seed(seed)
  fit <- abart(x.train = X, times = times, delta = delta, x.test = X_test, K = 2, ntree = 50L, k = 2, power = 2, base = 0.95,
               sigdf = 3, sigquant = 0.9, sigest = sigest, ndpost = ndpost, nskip = nskip, keepevery = 1L, printevery = 100000L)
  list(mu = fit$yhat.test, sigma = tail(fit$sigma, ndpost), sigest = if (is.na(sigest)) NA else sigest) }
one <- function(r) {
  out <- list()
  # repo DGP replicate
  d <- readRDS(file.path(SCR, "mapbart-sim-survival/data_v2", sprintf("data_p10_sc1E_alternative_%d.RData", r)))
  X <- as.matrix(d$X[, paste0("X", 1:10)]); D <- d$X$D; Z <- d$X$Z; ic <- Z == 0 & D == 1; it <- Z == 1; itest <- D == 1
  lower <- min(d$y[d$delta == 1]) / 10; eff <- exp(d$treat_eff_true)
  # our DGP replicate
  dat <- gen_data(1, "survival", seed = rep_seed(2026, r), n_mc = 0)
  Xo <- rbind(as.matrix(dat$rct_ctrl[, paste0("X", 1:10)])); Xt <- as.matrix(dat$rct_trt[, paste0("X", 1:10)])
  lower_o <- min(c(dat$rct_ctrl$time[dat$rct_ctrl$status == 1], dat$rct_trt$time[dat$rct_trt$status == 1], dat$rwd_ctrl$time[dat$rwd_ctrl$status == 1])) / 10
  for (mode in c("survreg", "default")) {
    fc <- fit_arm(X[ic, ], d$y[ic], d$delta[ic], X[itest, ], 1000 + r, mode); ft <- fit_arm(X[it, ], d$y[it], d$delta[it], X[itest, ], 2000 + r, mode)
    ps <- repo_pop_median(ft$mu, ft$sigma, lower) / repo_pop_median(fc$mu, fc$sigma, lower); ps <- ps[is.finite(ps)]
    out[[length(out) + 1]] <- data.frame(dgp = "repo", sigest = mode, rep = r, med = median(ps) - eff, mean = mean(ps) - eff, sd = sd(ps),
      cov = as.numeric(quantile(ps, .025) <= eff & quantile(ps, .975) >= eff), sig_c = mean(fc$sigma), sigest_c = fc$sigest)
    fc <- fit_arm(Xo, dat$rct_ctrl$time, dat$rct_ctrl$status, dat$x_rct, 1000 + r, mode); ft <- fit_arm(Xt, dat$rct_trt$time, dat$rct_trt$status, dat$x_rct, 2000 + r, mode)
    ps <- repo_pop_median(ft$mu, ft$sigma, lower_o) / repo_pop_median(fc$mu, fc$sigma, lower_o); ps <- ps[is.finite(ps)]; eff2 <- dat$truth$ate_med
    out[[length(out) + 1]] <- data.frame(dgp = "ours", sigest = mode, rep = r, med = median(ps) - eff2, mean = mean(ps) - eff2, sd = sd(ps),
      cov = as.numeric(quantile(ps, .025) <= eff2 & quantile(ps, .975) >= eff2), sig_c = mean(fc$sigma), sigest_c = fc$sigest)
  }
  do.call(rbind, out) }
res <- do.call(rbind, mclapply(1:R, one, mc.cores = 12))
ag <- do.call(rbind, lapply(split(res, list(res$dgp, res$sigest)), function(s) data.frame(dgp = s$dgp[1], sigest = s$sigest[1], R = nrow(s),
  bias_median = mean(s$med), bias_mean = mean(s$mean), se_bias = sd(s$med)/sqrt(nrow(s)), sd = mean(s$sd),
  rmse_median = sqrt(mean(s$med^2)), rmse_mean = sqrt(mean(s$mean^2)), cov = mean(s$cov), post_sigma_ctrl = mean(s$sig_c), sigest_ctrl = mean(s$sigest_c, na.rm = TRUE))))
print(ag, digits = 3, row.names = FALSE)
# paired effect of sigest on the same data
w <- reshape(res[, c("dgp", "sigest", "rep", "med")], idvar = c("dgp", "rep"), timevar = "sigest", direction = "wide")
cat("paired (survreg - default) median-bias: repo dgp", round(mean((w$med.survreg - w$med.default)[w$dgp == "repo"]), 3),
    " ours dgp", round(mean((w$med.survreg - w$med.default)[w$dgp == "ours"]), 3), "\n")
saveRDS(list(res = res, agg = ag), file.path(dirname(SCR), "recon/out/q3b_paired.rds"))
