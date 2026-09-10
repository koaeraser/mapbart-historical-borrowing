# check_hier_stan.R: cross-check the conjugate Gibbs sampler for HierLM
# (gibbs_hier in R/comparators.R) against the HMC fit of stan/hier_lm.stan on
# one Sc2 replicate. Both target the same posterior; posterior means and SDs
# of theta_1 (RCT control coefficients) and sigma_1 should agree to Monte
# Carlo error. Run from the 01-code root: Rscript tests/check_hier_stan.R
CODE_ROOT <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
suppressPackageStartupMessages(library(rstan))

dat <- gen_data_id("Sc2", "gaussian", seed = 5, n_mc = 0)
d1 <- prep_xy(dat$rct_ctrl, "gaussian"); d2 <- prep_xy(dat$rwd_ctrl, "gaussian")

set.seed(1)
gb <- gibbs_hier(d1$y, d1$X, d2$y, d2$X, n_burn = 2000, n_draw = 8000)

lambda <- lambda_bart(c(sigma_prelim(d1$y, d1$X), sigma_prelim(d2$y, d2$X)))
sd_ <- list(N = length(d1$y) + length(d2$y), K = ncol(d1$X) + 1,
            X = rbind(cbind(1, d1$X), cbind(1, d2$X)), y = c(d1$y, d2$y),
            g = c(rep(1L, length(d1$y)), rep(2L, length(d2$y))),
            lambda = lambda, nu = 3)
m <- stan_model(file.path(CODE_ROOT, "stan", "hier_lm.stan"))
fit <- sampling(m, data = sd_, chains = 4, iter = 4000, refresh = 0, seed = 1)
th <- extract(fit, "theta")$theta[, 1, ]
sg <- sqrt(extract(fit, "sigma2")$sigma2[, 1])

cmp <- data.frame(
  par = c(paste0("theta1_", 0:10), "sigma1"),
  gibbs_mean = c(colMeans(gb$theta), mean(gb$sigma)),
  stan_mean = c(colMeans(th), mean(sg)),
  gibbs_sd = c(apply(gb$theta, 2, sd), sd(gb$sigma)),
  stan_sd = c(apply(th, 2, sd), sd(sg)))
cmp$z <- (cmp$gibbs_mean - cmp$stan_mean) / cmp$stan_sd
print(cmp, digits = 3)
cat("max |mean diff| / posterior SD:", round(max(abs(cmp$z)), 3), "\n")
cat("max SD ratio deviation:", round(max(abs(cmp$gibbs_sd / cmp$stan_sd - 1)), 3), "\n")
if (max(abs(cmp$z)) < 0.15 && max(abs(cmp$gibbs_sd / cmp$stan_sd - 1)) < 0.15) {
  cat("check_hier_stan.R: PASS (Gibbs and HMC agree)\n")
} else {
  cat("check_hier_stan.R: FAIL\n")
}
