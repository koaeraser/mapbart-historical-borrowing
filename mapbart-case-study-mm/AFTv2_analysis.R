rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop"
projDir <- file.path(mainDir, "mapbart-case-study-mm")

library(rstan)
library(survival)
rstan_options(auto_write = TRUE)

OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
cat(sprintf("=== AFT_analysis.R: OUTCOME = %s ===\n", OUTCOME))

# ============================================================
# Single-layer lognormal AFT, mirroring mapbart-sim-survival/AFTv1.R:
#   - Control arm fit on UCMM (s=0) only
#   - Treatment arm fit on EloKRd (s=1) only
#   - No theta_1 / theta_2 hierarchy, no discrepancy variances
#   - BART-style priors (Chipman sigest -> lambda; total prior
#     variance range(y)^2/(4k^2), split w_alpha:w_beta = 2:1)
#   - mu_alpha = 0 (matching AFTv1)
#   - Estimand: RMST(5yr) ratio (EloKRd/UCMM) -- the restricted mean
#     survival time at tau = 5 yr, averaged over the EloKRd cohort,
#     under a uniform mixture of per-patient lognormals.  Matches the
#     HierAFT / mBART estimand (a bounded integral, far less
#     extrapolation-sensitive than a median ratio).
# ============================================================

# ============================================================
# 1. LOAD MERGED DATA
# ============================================================
# Select the merged cohort by FILE NAME via MERGED_FILE; default is the
# primary cohort (all regimens; E-Rd excluded). The n<N> token from the file name tags all outputs.
merged_file <- Sys.getenv("MERGED_FILE",
                          unset = file.path(projDir, "data_cleaned/merged_elokrd_ucmm_n283.RData"))
# Allow a bare filename: resolve against data_cleaned/.
if (!file.exists(merged_file) &&
    file.exists(file.path(projDir, "data_cleaned", basename(merged_file))))
  merged_file <- file.path(projDir, "data_cleaned", basename(merged_file))
if (!file.exists(merged_file))
  stop("Merged file not found: ", merged_file, "  (set MERGED_FILE; run data_merge.R first)")
load(merged_file)
data_tag <- regmatches(basename(merged_file), regexpr("n[0-9]+", basename(merged_file)))
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(merged_file))
cat(sprintf("Merged file: %s  (tag=%s)\n", basename(merged_file), data_tag))
cat(sprintf("Merged data: %d subjects (%d EloKRD, %d UCMM)\n",
            nrow(merged), sum(merged$trt == 1), sum(merged$trt == 0)))

# ============================================================
# 2. PREPARE COVARIATES
# ============================================================
covariate_cols <- c("age", "sex", "race", "ethnicity", "high_risk_cyto", "transplant_off_protocol")

prep_covariates <- function(df, cols) {
  X <- data.frame(row.names = 1:nrow(df))
  for (v in cols) {
    vals <- df[[v]]
    if (is.numeric(vals)) {
      X[[v]] <- vals
    } else {
      vals <- as.factor(vals)
      lvls <- levels(vals)
      if (length(lvls) == 2) {
        X[[v]] <- as.integer(vals == lvls[2])
      } else {
        for (k in 2:length(lvls)) {
          X[[paste0(v, "_", gsub("[/ ]", "_", lvls[k]))]] <- as.integer(vals == lvls[k])
        }
      }
    }
  }
  X
}

X_all <- prep_covariates(merged, covariate_cols)
cat(sprintf("Covariate matrix: %d x %d\n", nrow(X_all), ncol(X_all)))
cat("Columns:", paste(names(X_all), collapse = ", "), "\n")

trt <- merged$trt
if (OUTCOME == "PFS") {
  Y     <- merged$pfs_months
  event <- merged$pfs_status
} else {
  Y     <- merged$os_months
  event <- merged$os_status
}
# Convert months -> years so the log-survival scale matches the
# AFTv1 simulation (Y in years, mu_alpha = 0 reasonable).
Y <- Y / 12

trt_idx  <- which(trt == 1)
ctrl_idx <- which(trt == 0)

cat(sprintf("\nEloKRD (treatment): n=%d, events=%d, median %s=%.2f yr\n",
            length(trt_idx), sum(event[trt_idx]), OUTCOME, median(Y[trt_idx])))
cat(sprintf("UCMM (ext control): n=%d, events=%d, median %s=%.2f yr\n",
            length(ctrl_idx), sum(event[ctrl_idx]), OUTCOME, median(Y[ctrl_idx])))

# ============================================================
# 3. STAN MODEL (single-arm AFT, used for both arms)
# ============================================================
aft_stan_code <- "
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] y;
  int<lower=0,upper=1> event[N];
  vector<lower=0>[N] wt;       // per-obs likelihood weight (power-prior downweight of external control)
  real<lower=0> sigma_scale;
  real mu_alpha;
  real<lower=0> lambda_alpha;
  real<lower=0> lambda_beta;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> sigma2;
}
model {
  alpha ~ normal(mu_alpha, lambda_alpha);
  beta  ~ normal(0, lambda_beta);
  sigma2 ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0);
  {
    vector[N] mu = alpha + X * beta;
    real sigma = sqrt(sigma2);
    // power-prior weighted: each row enters the likelihood with weight wt[n]
    for (n in 1:N) {
      if (event[n] == 1) target += wt[n] * normal_lpdf(y[n] | mu[n], sigma);
      else               target += wt[n] * normal_lccdf(y[n] | mu[n], sigma);
    }
  }
}
"

cat("\nCompiling Stan AFT model...\n")
aft_model <- stan_model(model_code = aft_stan_code)

# ============================================================
# 4. MCMC + PRIOR SETTINGS (matched to AFTv1.R)
# ============================================================
n_chains <- 1
n_iter   <- 1600
n_warmup <- 100
n_cores  <- 3
set.seed(42)

nu       <- 3
sigquant <- 0.90
qchi     <- qchisq(1 - sigquant, nu)
k        <- 2.0
w_alpha  <- 2.0
w_beta   <- 1.0

# ---------------------------------------------------------------
# Power-prior weight on the EXTERNAL CONTROL (UCMM) likelihood.
# Mirrors mapbart-sim-survival-single-arm/AFTv2.R's rwd_w: each UCMM control row
# enters the control-arm likelihood with weight rwd_w, so the effective
# external-control sample size is ~ rwd_w * n_UCMM.  To match a MAP-BART
# borrowing target N: rwd_w = N / n_UCMM.  rwd_w = 1 reproduces full
# (unweighted) borrowing.  The treatment arm always gets weight 1.
# AFTV2_RWD_W = comma-separated weights to sweep (one result file per
# value, tagged _w<rwd_w>); unset -> full borrowing only.
.rwd_w_env <- strsplit(Sys.getenv("AFTV2_RWD_W", unset = "1"), "[, ]+")[[1]]
rwd_w_vals <- as.numeric(.rwd_w_env[nzchar(.rwd_w_env)])
rwd_w_vals <- rwd_w_vals[!is.na(rwd_w_vals)]
stopifnot(length(rwd_w_vals) >= 1)
n_rwd <- length(ctrl_idx)   # external-control (UCMM) count, for matched-N labeling
cat(sprintf("AFTv2 rwd_w sweep: %s  (n_UCMM = %d)\n",
            paste(rwd_w_vals, collapse = ", "), n_rwd))

# Center all covariates at the EloKRd (target population) means
x_means <- colMeans(as.matrix(X_all[trt_idx, ]))
## Do not center binary (0/1) covariates -- only continuous predictors are centered.
.is_binary <- apply(as.matrix(X_all), 2, function(col) all(col[!is.na(col)] %in% c(0, 1)))
x_means[.is_binary] <- 0

fit_one_arm <- function(idx, label, seed, wt = 1) {
  x <- as.matrix(X_all[idx, ])
  x <- scale(x, center = x_means, scale = FALSE)
  y <- log(Y[idx])
  e <- as.integer(event[idx])

  sds <- apply(x, 2, sd)
  x_aft <- x[, sds > 1e-8, drop = FALSE]
  if (ncol(x_aft) == 0) {
    fit_init <- survreg(Surv(exp(y), e) ~ 1, dist = "lognormal")
  } else {
    fit_init <- survreg(Surv(exp(y), e) ~ ., data = data.frame(x_aft),
                        dist = "lognormal")
  }
  sigest <- fit_init$scale
  lambda <- (sigest^2 * qchi) / nu

  total_var    <- (max(y) - min(y))^2 / (2 * k)^2
  total_weight <- w_alpha + ncol(x) * w_beta
  lambda_alpha    <- sqrt(total_var * w_alpha / total_weight)
  lambda_beta     <- sqrt(total_var * w_beta  / total_weight)

  stan_data <- list(
    N = length(y), P = ncol(x), X = x, y = y, event = e,
    wt = rep(wt, length(y)),   # power-prior weight (1 = full; <1 = downweighted)
    sigma_scale    = lambda,
    mu_alpha = 0,           # matches AFTv1
    lambda_alpha      = lambda_alpha,
    lambda_beta       = lambda_beta
  )

  cat(sprintf("  [%s] n=%d, p=%d, sigest=%.3f, lambda_alpha=%.3f, lambda_beta=%.3f, wt=%g\n",
              label, length(y), ncol(x), sigest, lambda_alpha, lambda_beta, wt))
  fit <- sampling(aft_model, data = stan_data,
                  chains = n_chains, iter = n_iter, warmup = n_warmup,
                  cores = n_cores, seed = seed, refresh = 0,
                  control = list(adapt_delta = 0.90))
  list(fit = fit, samples = rstan::extract(fit))
}

# ============================================================
# 5. FIT BOTH ARMS  (swept over the power-prior weight rwd_w)
# ============================================================
# Treatment arm is independent of rwd_w (always full weight), so fit it
# once; the control arm is refit per rwd_w with its likelihood downweighted.
cat("\n=== Fitting TREATMENT arm model (AFT on EloKRD) ===\n")
res_trt <- fit_one_arm(trt_idx, "treatment (EloKRd)", seed = 43L, wt = 1)
samples_trt <- res_trt$samples

for (rwd_w in rwd_w_vals) {
target_N <- as.integer(round(rwd_w * n_rwd))   # matched effective external-control N
cat(sprintf("\n########## AFTv2 rwd_w = %g  (target_N = %d / %d) ##########\n",
            rwd_w, target_N, n_rwd))

cat("\n=== Fitting CONTROL arm model (AFT on UCMM) ===\n")
res_ctrl <- fit_one_arm(ctrl_idx, "control (UCMM)", seed = 42L, wt = rwd_w)
samples_ctrl <- res_ctrl$samples

# ============================================================
# 6. PREDICT AT EloKRD COVARIATE PROFILE
# ============================================================
x_test <- as.matrix(X_all[trt_idx, ])
x_test <- scale(x_test, center = x_means, scale = FALSE)
n_test <- nrow(x_test)

n_samples_c <- length(samples_ctrl$alpha)
mu_pred_ctrl <- matrix(NA, n_samples_c, n_test)
for (i in 1:n_samples_c) {
  mu_pred_ctrl[i, ] <- samples_ctrl$alpha[i] + x_test %*% samples_ctrl$beta[i, ]
}
sig_draw_ctrl <- sqrt(samples_ctrl$sigma2)

n_samples_t <- length(samples_trt$alpha)
mu_pred_trt <- matrix(NA, n_samples_t, n_test)
for (i in 1:n_samples_t) {
  mu_pred_trt[i, ] <- samples_trt$alpha[i] + x_test %*% samples_trt$beta[i, ]
}
sig_draw_trt <- sqrt(samples_trt$sigma2)

cat(sprintf("\nControl  sigma: posterior mean=%.3f\n", mean(sig_draw_ctrl)))
cat(sprintf("Treatment sigma: posterior mean=%.3f\n", mean(sig_draw_trt)))

# ============================================================
# 7. POPULATION RMST (uniform mixture over EloKRd)
# ============================================================
# Primary estimand: RMST ratio at tau = 5 years -- matches HierAFT /
# mBART.  tau sits inside the observed follow-up so the integral is
# data-supported, avoiding the tail-extrapolation sensitivity of the
# median-survival estimand.
tau_rmst <- 5    # years (Y is in years)

# Closed form for E[min(Y, tau)] under lognormal(mu, sigma^2):
#   E[min(Y, tau)] = exp(mu + sigma^2/2) * Phi((log tau - mu - sigma^2)/sigma)
#                   + tau * Phi((mu - log tau)/sigma)
# Marginal RMST is the per-patient average over the EloKRd cohort.
compute_pop_rmst <- function(mu_draws, sig_draw, tau, label) {
  M  <- nrow(mu_draws)
  lt <- log(tau)
  rmst <- numeric(M)
  for (m in seq_len(M)) {
    mu <- mu_draws[m, ]; s <- sig_draw[m]
    rmst[m] <- mean(exp(mu + s^2 / 2) * pnorm((lt - mu - s^2) / s) +
                    tau * pnorm((mu - lt) / s))
  }
  cat(sprintf("  %s: posterior median RMST(%g yr) = %.2f yr\n",
              label, tau, median(rmst)))
  rmst
}

cat(sprintf("\n=== Computing population RMST at tau = %g years ===\n", tau_rmst))
rmst_ctrl <- compute_pop_rmst(mu_pred_ctrl, sig_draw_ctrl, tau_rmst, "Control (UCMM->EloKRD)")
rmst_trt  <- compute_pop_rmst(mu_pred_trt,  sig_draw_trt,  tau_rmst, "Treatment (EloKRD)")

# ============================================================
# 8. POSTERIOR TREATMENT EFFECT (RMST(5yr) ratio)
# ============================================================
cat(sprintf("\n=== Posterior Treatment Effect (RMST(%g yr) ratio) ===\n", tau_rmst))
post_ratio <- rmst_trt / rmst_ctrl
valid <- !is.na(post_ratio)
cat(sprintf("Valid posterior draws: %d / %d\n", sum(valid), length(post_ratio)))
post_ratio_valid <- post_ratio[valid]

delta_hat <- median(post_ratio_valid)
ci_95 <- quantile(post_ratio_valid, probs = c(0.025, 0.975))

cat(sprintf("\n  Posterior median RMST ratio (trt/ctrl): %.3f\n", delta_hat))
cat(sprintf("  95%% credible interval: [%.3f, %.3f]\n", ci_95[1], ci_95[2]))

# --- RMST difference (trt - ctrl): stable when control RMST -> 0 ---
post_diff_valid <- (rmst_trt - rmst_ctrl)[valid]
delta_diff <- median(post_diff_valid, na.rm = TRUE)
ci_diff_95 <- quantile(post_diff_valid, probs = c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("  Posterior median RMST difference (trt-ctrl): %.3f yr  [%.3f, %.3f]\n",
            delta_diff, ci_diff_95[1], ci_diff_95[2]))

# ============================================================
# 8b. ARM-SPECIFIC ESTIMATES (population RMST + residual sigma, median [95% CI])
# ============================================================
arm_summ <- function(x) {
  q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
  c(est = median(x, na.rm = TRUE), lo = unname(q[1]), hi = unname(q[2]))
}
rmst_trt_est   <- arm_summ(rmst_trt)
rmst_ctrl_est  <- arm_summ(rmst_ctrl)
sigma_trt_est  <- arm_summ(sig_draw_trt)
sigma_ctrl_est <- arm_summ(sig_draw_ctrl)
cat("\n=== Arm-specific estimates (median [95% CI]) ===\n")
cat(sprintf("  RMST  treatment: %.3f  [%.3f, %.3f]\n", rmst_trt_est["est"],  rmst_trt_est["lo"],  rmst_trt_est["hi"]))
cat(sprintf("  RMST  control  : %.3f  [%.3f, %.3f]\n", rmst_ctrl_est["est"], rmst_ctrl_est["lo"], rmst_ctrl_est["hi"]))
cat(sprintf("  sigma treatment: %.3f  [%.3f, %.3f]\n", sigma_trt_est["est"],  sigma_trt_est["lo"],  sigma_trt_est["hi"]))
cat(sprintf("  sigma control  : %.3f  [%.3f, %.3f]\n", sigma_ctrl_est["est"], sigma_ctrl_est["lo"], sigma_ctrl_est["hi"]))

# ============================================================
# 9. SAVE RESULTS
# ============================================================
out_dir <- file.path(projDir, "res")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
results <- list(
  post_ratio    = post_ratio_valid,
  delta_hat     = delta_hat,
  ci_95         = ci_95,
  delta_diff    = delta_diff,
  ci_diff_95    = ci_diff_95,
  tau_rmst      = tau_rmst,
  rmst_ctrl     = rmst_ctrl,
  rmst_trt      = rmst_trt,
  mu_pred_ctrl  = mu_pred_ctrl,
  mu_pred_trt   = mu_pred_trt,
  sig_draw_ctrl = sig_draw_ctrl,
  sig_draw_trt  = sig_draw_trt,
  rmst_trt_est   = rmst_trt_est,
  rmst_ctrl_est  = rmst_ctrl_est,
  sigma_trt_est  = sigma_trt_est,
  sigma_ctrl_est = sigma_ctrl_est,
  samples_ctrl  = samples_ctrl,
  samples_trt   = samples_trt,
  settings      = list(n_chains = n_chains, n_iter = n_iter,
                       n_warmup = n_warmup, nu = nu, sigquant = sigquant,
                       k = k, w_alpha = w_alpha, w_beta = w_beta,
                       rwd_w = rwd_w, n_rwd = n_rwd, target_N = target_N)
)
out_file <- file.path(out_dir, sprintf("AFTv2_results_%s_%s_w%g.RData",
                                       OUTCOME, data_tag, rwd_w))
saveRDS(results, file = out_file)
cat(sprintf("\nResults saved to: %s\n", out_file))

}  # end loop over rwd_w_vals

cat("\n=== DONE ===\n")
