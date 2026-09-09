rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop"
projDir <- file.path(mainDir, "mapbart-case-study-mm")

library(rstan)
library(survival)
rstan_options(auto_write = TRUE)   # cache compiled Stan models between runs

OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
cat(sprintf("=== HierAFT_analysis.R: OUTCOME = %s ===\n", OUTCOME))

# ============================================================
# AFT with (alpha, beta, sigma^2) hierarchy echoing
# MAP-BART's per-leaf (theta_1, theta_2, tau^2) structure.
#
# Control arm carries two parameter layers:
#   - External (theta_2 analog): (alpha2, beta2) fit on s=0 (UCMM).
#   - Concurrent (theta_1 analog): (alpha1, beta1) centered on the
#     external layer via alpha1 | alpha2 ~ N(alpha2, tau2_alpha),
#     beta1 | beta2 ~ N(beta2, tau2_beta); no s=1 data to fit.
#   - Shared residual variance sigma^2 across both layers, since the
#     s=1-specific residual variance is unidentified without concurrent
#     control data.
#
# Predictions for EloKRd test points use the concurrent layer.
# With no s=1 data, the concurrent layer is prior-informed only
# through the hierarchy, matching mBART's behavior when no
# concurrent-control data is available.
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

trt   <- merged$trt           # 1 = EloKRD, 0 = UCMM
if (OUTCOME == "PFS") {
  Y     <- merged$pfs_months
  event <- merged$pfs_status
} else {
  Y     <- merged$os_months
  event <- merged$os_status
}
# Convert months -> years so the log-survival scale matches the other analyses.
Y <- Y / 12

trt_idx  <- which(trt == 1)
ctrl_idx <- which(trt == 0)

cat(sprintf("\nEloKRD (treatment): n=%d, events=%d, median %s=%.1f mo\n",
            length(trt_idx), sum(event[trt_idx]), OUTCOME, median(Y[trt_idx])))
cat(sprintf("UCMM (ext control): n=%d, events=%d, median %s=%.1f mo\n",
            length(ctrl_idx), sum(event[ctrl_idx]), OUTCOME, median(Y[ctrl_idx])))

# ============================================================
# 3. STAN MODELS
# ============================================================

# Control arm: theta_1 / theta_2 hierarchy.
# External layer fit on s=0 (UCMM); concurrent layer centered on
# external layer with coefficient-level discrepancy variances
# (tau2_alpha, tau2_beta ~ IG(3/2, 3*prior/2), prior fixed).
aft_stan_ctrl <- "
data {
  int<lower=0> N2;                 // s=0 (external) sample size
  int<lower=0> N1;                 // s=1 (concurrent) sample size
  int<lower=0> P;
  matrix[N2, P] X2;
  vector[N2] y2;
  int<lower=0,upper=1> event2[N2];
  matrix[N1, P] X1;
  vector[N1] y1;
  int<lower=0,upper=1> event1[N1];
  real<lower=0> sigma_scale;
  real mu_alpha;
  real<lower=0> lambda_alpha;
  real<lower=0> lambda_beta;
  real<lower=0> prior;             // scale for tau2 ~ inv_gamma(3/2, 3*prior/2)
}
parameters {
  real alpha2;                     // theta_2 intercept
  vector[P] beta2;                 // theta_2 slopes
  real<lower=0> sigma_sq;          // shared residual variance across s=0 and s=1
  real alpha1;                     // theta_1 intercept
  vector[P] beta1;                 // theta_1 slopes
  real<lower=0> tau2_alpha;
  real<lower=0> tau2_beta;
}
model {
  // External-layer priors (match baseline AFT prior)
  alpha2 ~ normal(mu_alpha, lambda_alpha);
  beta2  ~ normal(0, lambda_beta);
  sigma_sq ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0);

  // Discrepancy-variance priors: IG(3/2, 3*prior/2) (matching HierLM/HierAFT)
  tau2_alpha ~ inv_gamma(3.0/2.0, 3.0*prior/2.0);
  tau2_beta  ~ inv_gamma(3.0/2.0, 3.0*prior/2.0);

  // Concurrent layer centered on external layer (shared sigma_sq)
  alpha1 ~ normal(alpha2, sqrt(tau2_alpha));
  beta1  ~ normal(beta2, sqrt(tau2_beta));

  // Likelihood for external data (s=0)
  {
    vector[N2] mu2 = alpha2 + X2 * beta2;
    real sd_sh = sqrt(sigma_sq);
    for (n in 1:N2) {
      if (event2[n] == 1) target += normal_lpdf(y2[n] | mu2[n], sd_sh);
      else                target += normal_lccdf(y2[n] | mu2[n], sd_sh);
    }
  }
  // Likelihood for concurrent data (s=1) -- empty in this analysis
  if (N1 > 0) {
    vector[N1] mu1 = alpha1 + X1 * beta1;
    real sd_sh = sqrt(sigma_sq);
    for (n in 1:N1) {
      if (event1[n] == 1) target += normal_lpdf(y1[n] | mu1[n], sd_sh);
      else                target += normal_lccdf(y1[n] | mu1[n], sd_sh);
    }
  }
}
"

# Treatment arm: single-group AFT (fit on EloKRD only) -- unchanged.
aft_stan_trt <- "
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] y;
  int<lower=0,upper=1> event[N];
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
    for (n in 1:N) {
      if (event[n] == 1) target += normal_lpdf(y[n] | mu[n], sigma);
      else               target += normal_lccdf(y[n] | mu[n], sigma);
    }
  }
}
"

cat("\nCompiling Stan AFT models...\n")
aft_model_ctrl <- stan_model(model_code = aft_stan_ctrl)
aft_model_trt  <- stan_model(model_code = aft_stan_trt)

# ============================================================
# 4. MCMC SETTINGS  (matched to mBART_analysis.R)
# ============================================================
n_chains <- 1
n_iter   <- 1600
n_warmup <- 100
n_cores  <- 3
set.seed(42)

# BART-matched prior scales
nu       <- 3
sigquant <- 0.90
qchi     <- qchisq(1 - sigquant, nu)
k        <- 2.0
w_alpha  <- 2.0
w_beta   <- 1.0
ntree    <- 50    # mirrors mBART's number of trees (for discrepancy calibration)

## --- prior sweep for the discrepancy-variance IG(3/2, 3*prior/2) ----------
## HIERAFT_CONFIGS = comma-separated prior values, e.g. "0.05,0.5".
## Default mirrors mapbart-sim-survival-single-arm/run_all.R's lmv4_prior_vals = c(0.05, 0.5).
hieraft_priors <- as.numeric(strsplit(Sys.getenv("HIERAFT_CONFIGS", unset = "0.05,0.5"),
                                      "[, ]+")[[1]])
hieraft_priors <- hieraft_priors[!is.na(hieraft_priors)]
cat(sprintf("HierAFT prior sweep: %s\n", paste(hieraft_priors, collapse = ", ")))

for (prior in hieraft_priors) {
cat(sprintf("\n############ HierAFT prior = %g ############\n", prior))

# ============================================================
# 5. CONTROL ARM: theta_2 fit on UCMM, theta_1 centered on it
# ============================================================
cat("\n=== Fitting CONTROL arm model (theta_1/theta_2 hierarchy) ===\n")

x_means <- colMeans(as.matrix(X_all[trt_idx, ]))
## Do not center binary (0/1) covariates -- only continuous predictors are centered.
.is_binary <- apply(as.matrix(X_all), 2, function(col) all(col[!is.na(col)] %in% c(0, 1)))
x_means[.is_binary] <- 0

x_ext <- as.matrix(X_all[ctrl_idx, ])
x_ext <- scale(x_ext, center = x_means, scale = FALSE)
y_ext <- log(Y[ctrl_idx])
event_ext <- as.integer(event[ctrl_idx])

# Frequentist sigma estimate (Chipman-style)
ext_sds <- apply(x_ext, 2, sd)
x_ext_aft <- x_ext[, ext_sds > 1e-8, drop = FALSE]
if (ncol(x_ext_aft) == 0) {
  aft_fit_ext <- survreg(Surv(exp(y_ext), event_ext) ~ 1, dist = "lognormal")
} else {
  aft_fit_ext <- survreg(Surv(exp(y_ext), event_ext) ~ .,
                         data = data.frame(x_ext_aft),
                         dist = "lognormal")
}
sigest_ext  <- aft_fit_ext$scale
lambda_ext  <- (sigest_ext^2 * qchi) / nu

# Baseline prior variances (external layer) -- split by w_alpha / w_beta
total_var_ext    <- (max(y_ext) - min(y_ext))^2 / (2 * k)^2
total_weight_ext <- w_alpha + ncol(x_ext) * w_beta
var_alpha_ext    <- total_var_ext * w_alpha / total_weight_ext
var_beta_ext     <- total_var_ext * w_beta  / total_weight_ext
lambda_alpha_ext    <- sqrt(var_alpha_ext)
lambda_beta_ext     <- sqrt(var_beta_ext)

# Discrepancy-variance prior: IG(3/2, 3*prior/2); `prior` is the sweep-loop variable.
# IG(3/2, 3*prior/2) has mean 3*prior.

cat(sprintf("  External prior: alpha ~ N(%.2f, %.3f^2), beta ~ N(0, %.3f^2)\n",
            mean(y_ext), lambda_alpha_ext, lambda_beta_ext))
cat(sprintf("  Discrepancy: tau2_alpha, tau2_beta ~ IG(3/2, 3*%.3g/2)  [E=%.4f]\n",
            prior, 3 * prior))

stan_data_ctrl <- list(
  N2 = length(y_ext),
  N1 = 0L,
  P  = ncol(x_ext),
  X2 = x_ext,
  y2 = y_ext,
  event2 = event_ext,
  X1 = matrix(0, 0, ncol(x_ext)),
  y1 = numeric(0),
  event1 = integer(0),
  sigma_scale    = lambda_ext,
  mu_alpha = mean(y_ext),
  lambda_alpha      = lambda_alpha_ext,
  lambda_beta       = lambda_beta_ext,
  prior = prior
)

cat(sprintf("  Training: N2=%d (s=0), N1=%d (s=1), P=%d\n",
            stan_data_ctrl$N2, stan_data_ctrl$N1, stan_data_ctrl$P))
cat("  Running Stan MCMC...\n")

fit_ctrl <- sampling(aft_model_ctrl,
                     data = stan_data_ctrl,
                     chains = n_chains,
                     iter = n_iter,
                     warmup = n_warmup,
                     cores = n_cores,
                     seed = 42L,
                     refresh = 0,
                     control = list(adapt_delta = 0.90))

samples_ctrl <- rstan::extract(fit_ctrl)

# Predict counterfactual concurrent-control outcomes for EloKRD patients.
# Uses the concurrent layer (alpha1, beta1) with shared sigma^2 -- the "theta_1 analog".
x_test <- as.matrix(X_all[trt_idx, ])
x_test <- scale(x_test, center = x_means, scale = FALSE)

n_samples <- length(samples_ctrl$alpha1)
n_test    <- nrow(x_test)

mu_pred_ctrl <- matrix(NA, n_samples, n_test)
for (i in 1:n_samples) {
  mu_pred_ctrl[i, ] <- samples_ctrl$alpha1[i] +
                       x_test %*% samples_ctrl$beta1[i, ]
}
sig_draw_ctrl <- sqrt(samples_ctrl$sigma_sq)   # shared sigma^2 across s=0 and s=1

cat(sprintf("  Control model done. sigma (shared): mean=%.3f\n",
            mean(sig_draw_ctrl)))
cat(sprintf("  Posterior of discrepancy variances:\n"))
cat(sprintf("    tau2_alpha: mean=%.3f (prior E=%.3f)\n",
            mean(samples_ctrl$tau2_alpha), 3 * prior))
cat(sprintf("    tau2_beta : mean=%.3f (prior E=%.3f)\n",
            mean(samples_ctrl$tau2_beta),  3 * prior))

# ============================================================
# 6. TREATMENT ARM: unchanged
# ============================================================
cat("\n=== Fitting TREATMENT arm model (AFT on EloKRD) ===\n")

x_trt <- as.matrix(X_all[trt_idx, ])
x_trt <- scale(x_trt, center = x_means, scale = FALSE)
y_trt <- log(Y[trt_idx])
event_trt <- as.integer(event[trt_idx])

trt_sds <- apply(x_trt, 2, sd)
x_trt_aft <- x_trt[, trt_sds > 1e-8, drop = FALSE]
if (ncol(x_trt_aft) == 0) {
  aft_fit_trt <- survreg(Surv(exp(y_trt), event_trt) ~ 1, dist = "lognormal")
} else {
  aft_fit_trt <- survreg(Surv(exp(y_trt), event_trt) ~ .,
                         data = data.frame(x_trt_aft),
                         dist = "lognormal")
}
sigest_trt <- aft_fit_trt$scale
lambda_trt <- (sigest_trt^2 * qchi) / nu

total_var_trt    <- (max(y_trt) - min(y_trt))^2 / (2 * k)^2
total_weight_trt <- w_alpha + ncol(x_trt) * w_beta
var_alpha_trt    <- total_var_trt * w_alpha / total_weight_trt
var_beta_trt     <- total_var_trt * w_beta  / total_weight_trt
lambda_alpha_trt    <- sqrt(var_alpha_trt)
lambda_beta_trt     <- sqrt(var_beta_trt)

stan_data_trt <- list(
  N = length(y_trt),
  P = ncol(x_trt),
  X = x_trt,
  y = y_trt,
  event = event_trt,
  sigma_scale    = lambda_trt,
  mu_alpha = mean(y_trt),
  lambda_alpha      = lambda_alpha_trt,
  lambda_beta       = lambda_beta_trt
)

cat(sprintf("  Training: n=%d, p=%d\n", length(y_trt), ncol(x_trt)))
cat("  Running Stan MCMC...\n")

fit_trt <- sampling(aft_model_trt,
                    data = stan_data_trt,
                    chains = n_chains,
                    iter = n_iter,
                    warmup = n_warmup,
                    cores = n_cores,
                    seed = 43L,
                    refresh = 0,
                    control = list(adapt_delta = 0.90))

samples_trt <- rstan::extract(fit_trt)

mu_pred_trt <- matrix(NA, length(samples_trt$alpha), n_test)
for (i in 1:length(samples_trt$alpha)) {
  mu_pred_trt[i, ] <- samples_trt$alpha[i] + x_test %*% samples_trt$beta[i, ]
}
sig_draw_trt <- sqrt(samples_trt$sigma2)

cat(sprintf("  Treatment model done. sigma: mean=%.3f\n", mean(sig_draw_trt)))

# ============================================================
# 7. POPULATION RMST (mixture over EloKRD)
# ============================================================
# Primary estimand: RMST ratio at tau = 5 years (matches KM/AFTv2/BARTv2/mBART).
# tau is chosen inside the observed follow-up horizon so the integral is
# data-supported, avoiding the tail-extrapolation sensitivity of
# the median-survival estimand.
tau_rmst <- 5    # years (Y is in years)
cat(sprintf("\n=== Computing population RMST at tau = %g years ===\n", tau_rmst))

# Closed form for E[min(Y, tau)] under lognormal(mu, sigma^2):
#   E[min(Y, tau)] = exp(mu + sigma^2/2) * Phi((log tau - mu - sigma^2)/sigma)
#                   + tau * Phi((mu - log tau)/sigma)
# Marginal RMST is the per-patient average.
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

rmst_ctrl <- compute_pop_rmst(mu_pred_ctrl, sig_draw_ctrl, tau_rmst,
                              "Control (UCMM->EloKRD via theta_1)")
rmst_trt  <- compute_pop_rmst(mu_pred_trt,  sig_draw_trt,  tau_rmst,
                              "Treatment (EloKRD)")

# ============================================================
# 8. POSTERIOR TREATMENT EFFECT
# ============================================================
cat(sprintf("\n=== Posterior Treatment Effect (RMST(%g) ratio) ===\n", tau_rmst))

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
  rmst_ctrl     = rmst_ctrl,
  rmst_trt      = rmst_trt,
  tau_rmst      = tau_rmst,
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
  calibration   = list(prior      = prior,
                       E_tau2     = 3 * prior,   # IG(3/2, 3*prior/2) mean
                       ntree_ref  = ntree),
  settings      = list(n_chains = n_chains, n_iter = n_iter,
                       n_warmup = n_warmup, nu = nu, sigquant = sigquant,
                       k = k, w_alpha = w_alpha, w_beta = w_beta)
)
out_file <- file.path(out_dir, sprintf("HierAFT_results_%s_%s_prior%g.RData", OUTCOME, data_tag, prior))
saveRDS(results, file = out_file)
cat(sprintf("\nResults saved to: %s\n", out_file))

}  # end prior sweep

cat("\n=== DONE ===\n")
