# run_reference.R: reproduce the paper's reference rows (Table rmst-results)
# with Yunxuan's own model code: Kaplan-Meier RMST ratio (km_analysis.R),
# complete-pooling log-normal AFT (AFTv2_analysis.R at rwd_w = 1), hierarchical
# AFT (HierAFT_analysis.R at prior 0.05 and 0.5; the paper does not say which
# fixed scale its row used, so both are run), and standard BART. The paper's
# BARTv2_analysis.R and mBART_analysis.R depend on aBART/cabart.cpp and
# mBART.real/cmbart.cpp, which are not in the repo, so the standard-BART row is
# fitted with lrcbart::bart_fit (plain AFT-BART, 50 trees, validated against
# BART::wbart in the package tests) and the MAP-AFT-BART rows are not rerun.
# Stan settings follow the scripts: 1 chain, 1600 iterations, 100 warm-up.
# Point estimates: Yunxuan's scripts report the posterior median; we store the
# draws and report mean, median and the 2.5/97.5 percentiles.
# Every model-based row is run under both covariate codings (see prep_data.R).
# One core.
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({library(rstan); library(survival); library(survRM2); library(lrcbart)})
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
APP <- file.path(Sys.getenv("LRC_REVISION_ROOT", "revision"), "04-application")
source(file.path(APP, "code/prep_data.R"))
TAU <- 5

# ---- Stan models, verbatim from AFTv2_analysis.R and HierAFT_analysis.R ----
aft_stan_code <- "
data {
  int<lower=0> N; int<lower=0> P; matrix[N, P] X; vector[N] y; int<lower=0,upper=1> event[N];
  vector<lower=0>[N] wt; real<lower=0> sigma_scale; real mu_alpha; real<lower=0> lambda_alpha; real<lower=0> lambda_beta;
}
parameters { real alpha; vector[P] beta; real<lower=0> sigma2; }
model {
  alpha ~ normal(mu_alpha, lambda_alpha); beta ~ normal(0, lambda_beta);
  sigma2 ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0);
  { vector[N] mu = alpha + X * beta; real sigma = sqrt(sigma2);
    for (n in 1:N) { if (event[n] == 1) target += wt[n] * normal_lpdf(y[n] | mu[n], sigma);
                     else target += wt[n] * normal_lccdf(y[n] | mu[n], sigma); } }
}"
hier_stan_code <- "
data {
  int<lower=0> N2; int<lower=0> N1; int<lower=0> P; matrix[N2, P] X2; vector[N2] y2; int<lower=0,upper=1> event2[N2];
  matrix[N1, P] X1; vector[N1] y1; int<lower=0,upper=1> event1[N1];
  real<lower=0> sigma_scale; real mu_alpha; real<lower=0> lambda_alpha; real<lower=0> lambda_beta; real<lower=0> prior;
}
parameters { real alpha2; vector[P] beta2; real<lower=0> sigma_sq; real alpha1; vector[P] beta1; real<lower=0> tau2_alpha; real<lower=0> tau2_beta; }
model {
  alpha2 ~ normal(mu_alpha, lambda_alpha); beta2 ~ normal(0, lambda_beta);
  sigma_sq ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0);
  tau2_alpha ~ inv_gamma(3.0/2.0, 3.0*prior/2.0); tau2_beta ~ inv_gamma(3.0/2.0, 3.0*prior/2.0);
  alpha1 ~ normal(alpha2, sqrt(tau2_alpha)); beta1 ~ normal(beta2, sqrt(tau2_beta));
  { vector[N2] mu2 = alpha2 + X2 * beta2; real sd_sh = sqrt(sigma_sq);
    for (n in 1:N2) { if (event2[n] == 1) target += normal_lpdf(y2[n] | mu2[n], sd_sh);
                      else target += normal_lccdf(y2[n] | mu2[n], sd_sh); } }
  if (N1 > 0) { vector[N1] mu1 = alpha1 + X1 * beta1; real sd_sh = sqrt(sigma_sq);
    for (n in 1:N1) { if (event1[n] == 1) target += normal_lpdf(y1[n] | mu1[n], sd_sh);
                      else target += normal_lccdf(y1[n] | mu1[n], sd_sh); } }
}"
cat("compiling Stan models\n")
aft_model <- stan_model(model_code = aft_stan_code)
hier_model <- stan_model(model_code = hier_stan_code)

nu <- 3; sigquant <- 0.90; qchi <- qchisq(1 - sigquant, nu); k <- 2; w_alpha <- 2; w_beta <- 1
n_iter <- 1600; n_warmup <- 100

rmst_mat <- function(mu, s, tau) { lt <- log(tau); s <- matrix(s, nrow(mu), ncol(mu))
  exp(mu + s^2 / 2) * pnorm((lt - mu - s^2) / s) + tau * pnorm((mu - lt) / s) }
pop_median <- function(mu, s) lrcbart:::lnorm_pop_median_mat(mu, s)
summ <- function(r) c(mean = mean(r), median = median(r), lo = unname(quantile(r, 0.025)), hi = unname(quantile(r, 0.975)))

prior_scales <- function(x, y, e) {
  sds <- apply(x, 2, sd); xa <- x[, sds > 1e-8, drop = FALSE]
  fi <- if (ncol(xa) == 0) survreg(Surv(exp(y), e) ~ 1, dist = "lognormal") else
    survreg(Surv(exp(y), e) ~ ., data = data.frame(xa), dist = "lognormal")
  lambda <- (fi$scale^2 * qchi) / nu
  total_var <- (max(y) - min(y))^2 / (2 * k)^2; tw <- w_alpha + ncol(x) * w_beta
  list(lambda = lambda, la = sqrt(total_var * w_alpha / tw), lb = sqrt(total_var * w_beta / tw))
}

fit_aft_arm <- function(x, y, e, seed, mu_alpha = 0) {
  ps <- prior_scales(x, y, e)
  sd <- list(N = length(y), P = ncol(x), X = x, y = y, event = e, wt = rep(1, length(y)),
             sigma_scale = ps$lambda, mu_alpha = mu_alpha, lambda_alpha = ps$la, lambda_beta = ps$lb)
  f <- sampling(aft_model, data = sd, chains = 1, iter = n_iter, warmup = n_warmup, cores = 1,
                seed = seed, refresh = 0, control = list(adapt_delta = 0.90))
  rstan::extract(f)
}
pred_lin <- function(s, alpha, beta, xt) { M <- length(alpha); mu <- matrix(NA, M, nrow(xt))
  for (i in 1:M) mu[i, ] <- as.vector(alpha[i] + xt %*% beta[i, ]); mu }

run_coding <- function(dat, coding) {
  X <- as.matrix(if (coding == "orig") dat$X_orig else dat$X_harm)
  ti <- dat$trt_idx; ci <- dat$ctrl_idx
  x_means <- colMeans(X[ti, ]); isb <- apply(X, 2, function(v) all(v %in% c(0, 1))); x_means[isb] <- 0
  Xc <- scale(X, center = x_means, scale = FALSE)
  y <- log(dat$time); e <- dat$status
  xt <- Xc[ti, , drop = FALSE]
  out <- list()
  # complete-pooling AFT (AFTv2, rwd_w = 1)
  st <- fit_aft_arm(Xc[ti, ], y[ti], e[ti], seed = 43L, mu_alpha = 0)
  sc <- fit_aft_arm(Xc[ci, ], y[ci], e[ci], seed = 42L, mu_alpha = 0)
  mu1 <- pred_lin(st, st$alpha, st$beta, xt); mu0 <- pred_lin(sc, sc$alpha, sc$beta, xt)
  s1 <- as.vector(sqrt(st$sigma2)); s0 <- as.vector(sqrt(sc$sigma2))
  out[["AFT-CP"]] <- list(rmst = rowMeans(rmst_mat(mu1, s1, TAU)) / rowMeans(rmst_mat(mu0, s0, TAU)),
                          med = pop_median(mu1, s1) / pop_median(mu0, s0))
  # hierarchical AFT (HierAFT_analysis.R), treated arm mu_alpha = mean(y_trt)
  st2 <- fit_aft_arm(Xc[ti, ], y[ti], e[ti], seed = 43L, mu_alpha = mean(y[ti]))
  mu1h <- pred_lin(st2, st2$alpha, st2$beta, xt); s1h <- as.vector(sqrt(st2$sigma2))
  ps <- prior_scales(Xc[ci, ], y[ci], e[ci])
  for (pr in c(0.05, 0.5)) {
    sd <- list(N2 = length(ci), N1 = 0L, P = ncol(X), X2 = Xc[ci, ], y2 = y[ci], event2 = e[ci],
               X1 = matrix(0, 0, ncol(X)), y1 = numeric(0), event1 = integer(0),
               sigma_scale = ps$lambda, mu_alpha = mean(y[ci]), lambda_alpha = ps$la, lambda_beta = ps$lb, prior = pr)
    f <- sampling(hier_model, data = sd, chains = 1, iter = n_iter, warmup = n_warmup, cores = 1,
                  seed = 42L, refresh = 0, control = list(adapt_delta = 0.90))
    sh <- rstan::extract(f)
    mu0h <- pred_lin(sh, sh$alpha1, sh$beta1, xt); s0h <- as.vector(sqrt(sh$sigma_sq))
    out[[sprintf("HierAFT-%g", pr)]] <- list(
      rmst = rowMeans(rmst_mat(mu1h, s1h, TAU)) / rowMeans(rmst_mat(mu0h, s0h, TAU)),
      med = pop_median(mu1h, s1h) / pop_median(mu0h, s0h),
      tau2 = c(alpha = mean(sh$tau2_alpha), beta = mean(sh$tau2_beta)))
  }
  # standard BART (plain AFT-BART on each arm, 50 trees, lrcbart::bart_fit)
  bt <- bart_fit(NULL, X[ti, ], status = e[ti], time = dat$time[ti], H_f = 50, n_burn = 2000, n_draw = 2000, seed = 43L)
  bc <- bart_fit(NULL, X[ci, ], status = e[ci], time = dat$time[ci], H_f = 50, n_burn = 2000, n_draw = 2000, seed = 42L)
  a <- lrc_ate(bt, bc, X[ti, ], outcome = "survival", tau = TAU)
  out[["BART"]] <- list(rmst = a$ate_rmst, med = a$ate_med)
  out
}

res <- list(); rows <- list()
for (oc in c("PFS", "OS")) {
  dat <- load_mm(oc)
  # Kaplan-Meier (km_analysis.R): survRM2::rmst2 at tau = 5
  tau <- min(TAU, max(dat$time[dat$trt_idx]), max(dat$time[dat$ctrl_idx]))
  km <- rmst2(dat$time, dat$status, dat$trt, tau = tau)
  r <- km$unadjusted.result["RMST (arm=1)/(arm=0)", ]
  rows[[length(rows) + 1]] <- data.frame(outcome = oc, coding = "none", method = "KM", estimand = "rmst",
                                         mean = NA, median = unname(r["Est."]), lo = unname(r["lower .95"]), hi = unname(r["upper .95"]))
  for (cd in c("orig", "harm")) {
    cat(oc, cd, "\n")
    o <- run_coding(dat, cd)
    res[[paste(oc, cd)]] <- o
    for (m in names(o)) for (es in c("rmst", "med")) {
      s <- summ(o[[m]][[es]])
      rows[[length(rows) + 1]] <- data.frame(outcome = oc, coding = cd, method = m, estimand = es,
                                             mean = s["mean"], median = s["median"], lo = s["lo"], hi = s["hi"])
    }
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(APP, "results/reference_summary.csv"), row.names = FALSE)
saveRDS(res, file.path(APP, "results/reference_draws.rds"))
print(tab, digits = 3)
cat("done\n")
