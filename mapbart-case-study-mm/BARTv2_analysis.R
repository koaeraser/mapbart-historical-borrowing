rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop"
projDir <- file.path(mainDir, "mapbart-case-study-mm")

library(Rcpp)
library(RcppEigen)
library(survival)

# ============================================================
# Standard BART comparator (no borrowing), adapted from
# ../mapbart-sim-survival/BART.R.
#
# Difference from the simulation BART.R: that script carries a study
# indicator D (concurrent RCT control vs historical RWD control) as a
# covariate and trains the control fit on all Z==0 (D in {0,1}).  Here
# there is NO concurrent control and NO study indicator -- the control
# arm is purely external (UCMM).  So:
#   - Control BART: train on UCMM (trt==0), predict the EloKRd
#     counterfactual control outcome (test = EloKRd covariates).  No D.
#   - Treatment BART: train on EloKRd (trt==1), predict on EloKRd.
# Both arms use plain aBART (cabart) -- no theta_1/theta_2 hierarchy,
# no discrepancy variance, no calibrated s^2 (that is mBART).
#
# Estimand: RMST(5yr) ratio (EloKRd/UCMM), matching AFTv2/HierAFT/mBART/KM.
# Saves res/BARTv2_results_<OUTCOME>_<tag>.RData.
# ============================================================
OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
cat(sprintf("=== BARTv2_analysis.R: OUTCOME = %s ===\n", OUTCOME))

# ============================================================
# 1. LOAD MERGED DATA
# ============================================================
merged_file <- Sys.getenv("MERGED_FILE",
                          unset = file.path(projDir, "data_cleaned/merged_elokrd_ucmm_n283.RData"))
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
# 2. PREPARE DATA (covariates identical to mBART_analysis.R)
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

trt   <- merged$trt          # 1 = EloKRd, 0 = UCMM
if (OUTCOME == "PFS") {
  Y     <- merged$pfs_months
  event <- merged$pfs_status
} else {
  Y     <- merged$os_months
  event <- merged$os_status
}
Y <- Y / 12   # months -> years (mirror the other analyses; tau=5 yr)

trt_idx  <- which(trt == 1)   # EloKRd treatment arm
ctrl_idx <- which(trt == 0)   # UCMM external control

cat(sprintf("\nEloKRD (treatment): n=%d, events=%d, median %s=%.2f yr\n",
            length(trt_idx), sum(event[trt_idx]), OUTCOME, median(Y[trt_idx])))
cat(sprintf("UCMM (ext control): n=%d, events=%d, median %s=%.2f yr\n",
            length(ctrl_idx), sum(event[ctrl_idx]), OUTCOME, median(Y[ctrl_idx])))

# ============================================================
# 3. COMPILE aBART
# ============================================================
cat("\nCompiling aBART...\n")
sourceCpp(file.path(mainDir, "aBART/cabart.cpp"), rebuild = TRUE)
source(file.path(mainDir, "bartModelMatrix.R"))

# ============================================================
# 4. BART SETTINGS (default BART, both arms -- no shrinkage/borrowing)
# ============================================================
ndpost    <- 1500L
nskip     <- 200L
keepevery <- 1L
ntree     <- 50L
k         <- 2
numcut    <- 100L
alpha     <- 0.95   # base
beta_tree <- 2.0    # power
sigdf     <- 3
sigquant  <- 0.90

# ============================================================
# 5. FIT ONE BART ARM
#    train on `train_idx`, predict at `test_idx` covariates.
#    Returns mu draws (ndpost x n_test, log-survival scale) and sigma.
# ============================================================
fit_bart_arm <- function(train_idx, test_idx, seed, label) {
  cat(sprintf("\n=== Fitting %s (aBART) ===\n", label))
  x.train <- data.frame(X = as.matrix(X_all[train_idx, ]))
  y.train <- log(Y[train_idx])
  ev      <- as.integer(event[train_idx])
  x.test  <- data.frame(X = as.matrix(X_all[test_idx, ]))

  temp <- bartModelMatrix(x.train, numcut, usequants = FALSE,
                          xinfo = matrix(0, 0, 0), rm.const = TRUE)
  x.train.m <- t(temp$X)
  numcut.v  <- temp$numcut
  xinfo.m   <- temp$xinfo
  x.test.m  <- bartModelMatrix(x.test)
  x.test.m  <- t(x.test.m[, temp$rm.const])

  p  <- nrow(x.train.m); n <- ncol(x.train.m); np <- ncol(x.test.m)
  cat(sprintf("  Training: n=%d, p=%d; Test: np=%d\n", n, p, np))

  # Center on observed mean so BART learns residuals; cabart adds the
  # offset back so predictions stay on the original log-survival scale.
  offset <- mean(y.train)
  # offset <- 0   # uncentered version (commented out)
  y.sc   <- y.train - offset

  df  <- data.frame(t(x.train.m))
  aft <- survreg(Surv(exp(y.sc), ev) ~ ., data = df, dist = "lognormal")
  sigest <- aft$scale
  qchi   <- qchisq(1.0 - sigquant, sigdf)
  lambda <- (sigest^2 * qchi) / sigdf
  tau    <- (max(y.sc) - min(y.sc)) / (2 * k * sqrt(ntree))

  grp <- temp$grp; if (is.null(grp)) grp <- 1:p
  w   <- rep(1, n)

  res <- cabart(1L,            # ntype = wbart (continuous)
                n, p, np,
                x.train.m, y.sc, ev, x.test.m,
                ntree, numcut.v,
                ndpost * keepevery, nskip, keepevery,
                beta_tree,      # power
                alpha,          # base
                offset,
                tau, sigdf, lambda, sigest,
                w,
                FALSE,          # sparse
                0, 1,           # theta, omega
                grp,
                0.5, 1, p,      # a, b, rho
                FALSE,          # augment
                10000L,         # printevery
                xinfo.m,
                seed)

  # cabart's sigma includes the nskip burn-in draws; trim them (yhat.test
  # is already post-burn-in).  Mirrors mBART_analysis.R.
  if (nskip > 0) { nskip. <- 1:nskip } else { nskip. <- 0 }
  if (keepevery > 1)
    res$sigma <- c(res$sigma[nskip.],
                   res$sigma[nskip + seq(1, ndpost * keepevery, keepevery)])
  res$sigma <- res$sigma[-(nskip.)]

  cat(sprintf("  Done. sigma: mean=%.3f\n", mean(res$sigma)))
  list(mu = res$yhat.test, sigma = res$sigma)
}

set.seed(42)
# Control arm: train on UCMM, predict the EloKRd counterfactual control.
res_ctrl <- fit_bart_arm(ctrl_idx, trt_idx, seed = 42L,
                         "CONTROL arm (UCMM -> EloKRd counterfactual)")
# Treatment arm: train on EloKRd, predict on EloKRd.
res_trt  <- fit_bart_arm(trt_idx,  trt_idx, seed = 43L,
                         "TREATMENT arm (EloKRd)")

# ============================================================
# 6. POPULATION RMST (uniform mixture over EloKRd)
# ============================================================
tau_rmst <- 5    # years

# Closed form for E[min(Y, tau)] under lognormal(mu, sigma^2):
#   E[min(Y, tau)] = exp(mu + sigma^2/2) * Phi((log tau - mu - sigma^2)/sigma)
#                   + tau * Phi((mu - log tau)/sigma)
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
rmst_ctrl <- compute_pop_rmst(res_ctrl$mu, res_ctrl$sigma, tau_rmst, "Control (UCMM->EloKRD)")
rmst_trt  <- compute_pop_rmst(res_trt$mu,  res_trt$sigma,  tau_rmst, "Treatment (EloKRD)")

# ============================================================
# 7. POSTERIOR TREATMENT EFFECT (RMST(5yr) ratio)
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
# 7b. ARM-SPECIFIC ESTIMATES (population RMST + residual sigma, median [95% CI])
# ============================================================
arm_summ <- function(x) {
  q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
  c(est = median(x, na.rm = TRUE), lo = unname(q[1]), hi = unname(q[2]))
}
rmst_trt_est   <- arm_summ(rmst_trt)
rmst_ctrl_est  <- arm_summ(rmst_ctrl)
sigma_trt_est  <- arm_summ(res_trt$sigma)
sigma_ctrl_est <- arm_summ(res_ctrl$sigma)
cat("\n=== Arm-specific estimates (median [95% CI]) ===\n")
cat(sprintf("  RMST  treatment: %.3f  [%.3f, %.3f]\n", rmst_trt_est["est"],  rmst_trt_est["lo"],  rmst_trt_est["hi"]))
cat(sprintf("  RMST  control  : %.3f  [%.3f, %.3f]\n", rmst_ctrl_est["est"], rmst_ctrl_est["lo"], rmst_ctrl_est["hi"]))
cat(sprintf("  sigma treatment: %.3f  [%.3f, %.3f]\n", sigma_trt_est["est"],  sigma_trt_est["lo"],  sigma_trt_est["hi"]))
cat(sprintf("  sigma control  : %.3f  [%.3f, %.3f]\n", sigma_ctrl_est["est"], sigma_ctrl_est["lo"], sigma_ctrl_est["hi"]))

# ============================================================
# 8. SAVE RESULTS
# ============================================================
out_dir <- file.path(projDir, "res")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
results <- list(
  post_ratio = post_ratio_valid,
  delta_hat  = delta_hat,
  ci_95      = ci_95,
  delta_diff = delta_diff,
  ci_diff_95 = ci_diff_95,
  tau_rmst   = tau_rmst,
  rmst_ctrl  = rmst_ctrl,
  rmst_trt   = rmst_trt,
  rmst_trt_est   = rmst_trt_est,
  rmst_ctrl_est  = rmst_ctrl_est,
  sigma_trt_est  = sigma_trt_est,
  sigma_ctrl_est = sigma_ctrl_est,
  res_ctrl   = res_ctrl,
  res_trt    = res_trt,
  settings   = list(ndpost = ndpost, nskip = nskip, ntree = ntree, k = k,
                    alpha = alpha, beta = beta_tree, target_N = NA_integer_)
)
out_file <- file.path(out_dir, sprintf("BARTv2_results_%s_%s.RData", OUTCOME, data_tag))
saveRDS(results, file = out_file)
cat(sprintf("\nResults saved to: %s\n", out_file))

cat("\n=== DONE ===\n")
