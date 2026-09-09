rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop"
projDir <- file.path(mainDir, "mapbart-case-study-mm")

library(Rcpp)
library(RcppEigen)
library(survival)

# ============================================================
# OUTCOME FLAG ("PFS" default, or "OS")
#   PFS -> pfs_months/pfs_status;  OS -> os_months/os_status.
#   Results saved to res/mBART_results_<OUTCOME>_<tag>_N<target>.RData.
# The control-arm BART and the discrepancy-prior rate are the SAME for
# both outcomes: the shrunk control BART (ntree=5, k=10, see MCMC
# settings below), and the IG tau^2 rate = 3*prior_val/2 set at runtime
# from the calibrated s^2 (cmbart's setrate) -- no info.h edit/recompile.
# ============================================================
OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
cat(sprintf("=== mBART_analysis.R: OUTCOME = %s ===\n", OUTCOME))

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
# 2. PREPARE DATA FOR mBART
# ============================================================
# Covariates: all non-ID, non-outcome, non-date columns
covariate_cols <- c("age", "sex", "race", "ethnicity", "high_risk_cyto", "transplant_off_protocol")

# Convert categorical variables to numeric for BART
prep_covariates <- function(df, cols) {
  X <- data.frame(row.names = 1:nrow(df))
  for (v in cols) {
    vals <- df[[v]]
    if (is.numeric(vals)) {
      X[[v]] <- vals
    } else {
      # Convert to dummy variables
      vals <- as.factor(vals)
      lvls <- levels(vals)
      if (length(lvls) == 2) {
        # Binary: single indicator
        X[[v]] <- as.integer(vals == lvls[2])
      } else {
        # Multi-level: K-1 dummies
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

# Treatment and outcome
trt   <- merged$trt           # 1 = EloKRD, 0 = UCMM
if (OUTCOME == "PFS") {
  Y     <- merged$pfs_months
  event <- merged$pfs_status
} else {
  Y     <- merged$os_months
  event <- merged$os_status
}
# Convert months -> years to mirror AFT_analysis.R scaling.  No
# substantive effect (mBART centers Y via offset0 / trt_offset),
# but keeps internal log-time units consistent across analyses.
Y <- Y / 12

# Indices
trt_idx  <- which(trt == 1)   # EloKRD treatment arm
ctrl_idx <- which(trt == 0)   # UCMM external control

cat(sprintf("\nEloKRD (treatment): n=%d, events=%d, median %s=%.2f yr\n",
            length(trt_idx), sum(event[trt_idx]), OUTCOME,
            median(Y[trt_idx])))
cat(sprintf("UCMM (ext control): n=%d, events=%d, median %s=%.2f yr\n",
            length(ctrl_idx), sum(event[ctrl_idx]), OUTCOME,
            median(Y[ctrl_idx])))

# ============================================================
# 3. COMPILE mBART.real (modified for external-control-only)
# ============================================================
cat("\nCompiling mBART.real...\n")
sourceCpp(file.path(mainDir, "mBART.real/cmbart.cpp"), rebuild = TRUE)
cat("Compiling aBART...\n")
sourceCpp(file.path(mainDir, "aBART/cabart.cpp"), rebuild = TRUE)

# ============================================================
# 4. MCMC SETTINGS
# ============================================================
ndpost    <- 1500L   # matched to AFT_analysis.R for paired PP checks
nskip     <- 200L
keepevery <- 1L
# Control-arm BART hyperparameters.  The discrepancy-prior rate
# (3*prior_val/2) is set at runtime from the calibrated s^2, not
# compiled into info.h.
#
# Control-arm sensitivity sweep: the configs are NOT hard-coded here --
# they come from the MBART_CONFIGS env var (set by run_all.R, or
# directly), a comma-separated list of ntree:k pairs, e.g.
#   MBART_CONFIGS="5:10,20:5,50:2"
# Unset -> the full default sweep below.  Each config writes its own
# result files (tagged _nt<ntree>_k<k>) and adds a row to the comparison
# table printed at the end.
.cfg_spec <- Sys.getenv("MBART_CONFIGS", unset = "5:2,10:2,5:4")
control_configs <- lapply(strsplit(.cfg_spec, "[, ]+")[[1]], function(p) {
  if (!nzchar(p)) return(NULL)
  v <- strsplit(p, ":")[[1]]
  if (length(v) != 2L || anyNA(suppressWarnings(as.numeric(v))))
    stop("Bad MBART_CONFIGS entry '", p, "' (expected ntree:k, e.g. 5:10)")
  list(ntree = as.integer(v[1]), k = as.numeric(v[2]))
})
control_configs <- Filter(Negate(is.null), control_configs)
stopifnot(length(control_configs) >= 1)
cat(sprintf("Control-arm configs (MBART_CONFIGS): %s\n",
            paste(sapply(control_configs,
                         function(cf) sprintf("(ntree=%d,k=%g)", cf$ntree, cf$k)),
                  collapse = ", ")))
numcut    <- 100L
alpha     <- 0.95
beta_tree <- 2.0
sigdf     <- 3
sigquant  <- 0.90
# IG prior on per-leaf theta_1-theta_2 discrepancy variance tau^2.
# cmbart sets pi.rate = 3 * prior_val / 2 and pi.shape = 3/2, so
#   E[tau^2_discrepancy] = pi.rate / (pi.shape - 1) = 3 * prior_val.
#
# Build parallel vectors `prior_vals` and `target_Ns` from
# mapbart-case-study-mm/ess_calibration_surv.R's saved .RData -- one entry per
# N_target multiplier (0.50, 0.75, 1.00, 2.00, 3.00).  The analysis below
# loops over them and writes one result file per multiplier:
#   res/mBART_results_<OUTCOME>_N<target>.RData
# Multipliers > 1 (2x, 3x N_target) tighten the discrepancy (smaller s^2,
# more borrowing) -- used to probe the control-arm collapse.
#
# PRIOR_VAL_OVERRIDE=<num> env var bypasses the calibration entirely
# and runs a single iteration with that scalar (target_N := N_target
# from the .RData if available, else 100).
# Candidate pool; run_all.R extracts the ess-supported subset and passes it via
# MBART_NMULT.  Unset (standalone) -> the full pool, filtered by the interior-
# crossing rule in pick_prior_for below.
.nmult_env <- strsplit(Sys.getenv("MBART_NMULT", unset = "3,2,1,0.75,0.5"), "[, ]+")[[1]]
N_target_multipliers <- as.numeric(.nmult_env[nzchar(.nmult_env)])
.override <- suppressWarnings(as.numeric(Sys.getenv("PRIOR_VAL_OVERRIDE", "")))
# Fixed-s^2 MAP-BART (s^2 = 1e-4): a very tight discrepancy prior -> near-complete
# borrowing (the pooling extreme), independent of the ESS calibration.  Appended
# to every config's sweep and tagged N0.  MBART_FIXED_S2="" disables it.
.fixed_s2 <- suppressWarnings(as.numeric(Sys.getenv("MBART_FIXED_S2", "1e-4")))
.cal_dir  <- file.path(projDir, "ess_local", "res")

# Per-config calibration file: REQUIRE the exact config tag _nt<ntree>_k<k> so
# each control config reads its OWN ess calibration -- s^2_star depends on
# tau_theta = (max-min)/(2*k*sqrt(ntree)) and the leaf partial-residual structure,
# so a config's s^2 is not transferable to another.  Latest by mtime among files
# for that exact config, or NA.  Mirrors mapbart-sim-survival-realcov/mBART.R pick_prior_for.
cal_file_for <- function(ctrl_ntree, ctrl_k) {
  pat <- sprintf("^ess_calibration_ucmm_%s_.*_%s_nt%d_k%g\\.RData$",
                 tolower(OUTCOME), data_tag, ctrl_ntree, ctrl_k)
  fs  <- list.files(.cal_dir, pattern = pat, full.names = TRUE)
  fs  <- fs[!grepl("_checkpoint\\.RData$", fs)]
  if (!length(fs)) return(NA_character_)
  fs[order(file.info(fs)$mtime, decreasing = TRUE)][1]
}

# Build (prior_vals, target_Ns) for ONE control config from its config-matched
# calibration: same three branches (override / calibration / fallback) plus the
# fixed-s^2 (N0) append.  Called per config inside the loop below.
pick_prior_for <- function(ctrl_ntree, ctrl_k) {
  cf  <- cal_file_for(ctrl_ntree, ctrl_k)
  cal <- if (!is.na(cf)) readRDS(cf) else NULL
  if (is.finite(.override)) {
    pv <- .override
    tn <- if (!is.null(cal)) as.integer(round(cal$N_target)) else 100L
    cat(sprintf("  PRIOR_VAL_OVERRIDE -> prior_val = %.4g, target_N = %d\n", pv, tn))
  } else if (!is.null(cal)) {
    # Keep m only if Ghat(m x N_target) crosses q at an INTERIOR grid point.
    .pick <- function(m) {
      Ghat_m <- colMeans(cal$ESSbar <= m * cal$N_target, na.rm = TRUE)
      if (Ghat_m[1] >= cal$q || !any(Ghat_m >= cal$q)) return(NA_real_)
      cal$S_grid[min(which(Ghat_m >= cal$q))]
    }
    pv <- sapply(N_target_multipliers, .pick)
    tn <- as.integer(round(N_target_multipliers * cal$N_target))
    ok <- which(!is.na(pv)); pv <- pv[ok]; tn <- tn[ok]
    cat(sprintf("  calibration = %s\n", basename(cf)))
    for (j in seq_along(pv))
      cat(sprintf("    N=%d (%.2fxN) -> prior_val = %.4g\n",
                  tn[j], N_target_multipliers[ok][j], pv[j]))
  } else {
    pv <- 0.25; tn <- 100L
    cat(sprintf("  prior_val = %.4g, target_N = %d  (fallback; no %s/%s config-matched calibration)\n",
                pv, tn, tolower(OUTCOME), data_tag))
  }
  if (is.finite(.fixed_s2)) {
    pv <- c(pv, .fixed_s2); tn <- c(tn, 0L)
    cat(sprintf("  + fixed-s^2 MAP-BART -> prior_val (s^2) = %.4g (tagged N0)\n", .fixed_s2))
  }
  stopifnot(length(pv) == length(tn), length(pv) >= 1)
  list(prior_vals = pv, target_Ns = tn)
}

# Run only configs that HAVE a config-matched calibration (unless
# PRIOR_VAL_OVERRIDE bypasses calibration).  Mirrors mapbart-sim-survival-realcov/mBART.R.
if (!is.finite(.override)) {
  .has_cal <- vapply(control_configs, function(cf) !is.na(cal_file_for(cf$ntree, cf$k)), logical(1))
  if (any(!.has_cal))
    cat(sprintf("mBART_analysis.R: skipping configs with no %s/%s calibration: %s\n",
                tolower(OUTCOME), data_tag,
                paste(sapply(control_configs[!.has_cal],
                             function(cf) sprintf("(ntree=%d,k=%g)", cf$ntree, cf$k)), collapse = ", ")))
  control_configs <- control_configs[.has_cal]
  if (!length(control_configs))
    # No config calibrated -> skip MAP-BART entirely (loop is a no-op below),
    # rather than hard-stopping.  Mirrors mapbart-sim-survival-realcov/mBART.R.
    warning(sprintf("mBART_analysis.R: no control config has a matching %s/%s ess calibration (ess_calibration_ucmm_%s_*_%s_nt<ntree>_k<k>.RData in %s) -- SKIPPING MAP-BART. Run ess_calibration for the sweep's configs first.",
                    tolower(OUTCOME), data_tag, tolower(OUTCOME), data_tag, .cal_dir))
}

# =====================================================================
# Outer loop: control-arm BART config (ntree, k).  Inner loop: the
# N_target multipliers built above.  Each (config x target) runs the
# full analysis (control arm, treatment arm, posterior) and writes a
# result file tagged with the config and target_N.
# =====================================================================
mbart_summary <- list()   # one row per (config, target) for the final table

for (.cfg in control_configs) {
ntree <- .cfg$ntree
k     <- .cfg$k
cat(sprintf("\n############ CONTROL-ARM CONFIG: ntree=%d, k=%g ############\n",
            ntree, k))

# This config's own ess-selected s^2 sweep (config-matched; see pick_prior_for).
.pp        <- pick_prior_for(ntree, k)
prior_vals <- .pp$prior_vals
target_Ns  <- .pp$target_Ns

for (.idx in seq_along(prior_vals)) {

prior_val <- prior_vals[.idx]
target_N  <- target_Ns[.idx]
cat(sprintf("\n##### iteration %d/%d: N_target=%d, prior_val=%.4g #####\n",
            .idx, length(prior_vals), target_N, prior_val))

set.seed(42)

# ============================================================
# 5. CONTROL ARM: mBART on external control (UCMM only)
#    s=0 for all training data (external control)
#    s=0 for test data (predict for EloKRD patients under control)
#
#    With the modified bd.h, the model runs with only s=0 data,
#    using theta_2 for all predictions.
# ============================================================
cat("\n=== Fitting CONTROL arm model (mBART on UCMM) ===\n")

# Training: UCMM external control only
x.train.ctrl <- data.frame(X = as.matrix(X_all[ctrl_idx, ]))
s.train.ctrl <- rep(0, length(ctrl_idx))   # all external control (s=0)
y.train.ctrl <- log(Y[ctrl_idx])
event.ctrl   <- as.integer(event[ctrl_idx])

# Test: EloKRD patients (predict their counterfactual control outcome)
x.test.ctrl <- data.frame(X = as.matrix(X_all[trt_idx, ]))
s.test.ctrl <- rep(1, nrow(x.test.ctrl))   # predict using s=1 (concurrent control via theta_1)

# Prepare model matrix
source(file.path(mainDir, "bartModelMatrix.R"))

xinfo <- matrix(0, 0, 0)
temp <- bartModelMatrix(x.train.ctrl, numcut, usequants = FALSE,
                        xinfo = xinfo, rm.const = TRUE)
x.train.m <- t(temp$X)
numcut.v  <- temp$numcut
xinfo.m   <- temp$xinfo

x.test.m <- bartModelMatrix(x.test.ctrl)
x.test.m <- t(x.test.m[, temp$rm.const])

p_ctrl  <- nrow(x.train.m)
n_ctrl  <- ncol(x.train.m)
np_ctrl <- ncol(x.test.m)

cat(sprintf("  Training: n=%d, p=%d; Test: np=%d\n", n_ctrl, p_ctrl, np_ctrl))

# Sigma estimation via AFT.
# Center the s=0 data on its mean before fitting BART so the
# leaf prior N(0, lambda^2) is centered on the observed UCMM
# log-survival distribution rather than on zero.  The C++
# output layer adds offset0 back via `offset0 + bm.f * offset1`,
# so predictions remain on the original log-survival scale.
offset0 <- mean(y.train.ctrl)
# offset0 <- 0   # uncentered version (commented out)
offset1 <- 1
y.scaled <- (y.train.ctrl - offset0) / offset1

# Only external control data (s=0), so fit single AFT
df0 <- data.frame(t(x.train.m))
aft_fit0 <- survreg(Surv(exp(y.scaled), event.ctrl) ~ .,
                    data = df0, dist = "lognormal")
sigest0 <- aft_fit0$scale

# For the empty s=1 group, use same estimate as fallback
sigest1 <- sigest0

sigscale1 <- (sigest1^2) * qchisq(1 - sigquant, df = sigdf) / sigdf
sigscale0 <- (sigest0^2) * qchisq(1 - sigquant, df = sigdf) / sigdf

# ntree and k come from the control_configs loop (see MCMC settings).
tau <- (max(y.scaled) - min(y.scaled)) / (2 * k * sqrt(ntree))
cat(sprintf("  lambda = %.4f, tau_discrepancy prior_val = %.4f (E[tau^2] = %.4f)\n",
            tau, prior_val, 3 * prior_val))

type <- 1  # wbart (continuous)

cat("  Running mBART MCMC...\n")
res_ctrl <- cmbart(type,
                   n_ctrl,
                   p_ctrl,
                   np_ctrl,
                   x.train.m,
                   y.scaled,
                   s.train.ctrl,
                   event.ctrl,
                   x.test.m,
                   s.test.ctrl,
                   ntree,
                   numcut.v,
                   ndpost * keepevery,
                   nskip,
                   keepevery,
                   beta_tree,
                   alpha,
                   offset0,
                   offset1,
                   sigdf,
                   sigscale1,
                   sigscale0,
                   sigest1,
                   sigest0,
                   tau,
                   xinfo.m,
                   "IG",
                   prior_val,
                   42L)

# Extract sigma draws (post burn-in).  cmbart outputs a single sigma,
# sigma_rwd, shared across s=0 (UCMM) and s=1 (concurrent) groups,
# since with no concurrent-control data the s=1-specific residual
# variance is unidentified.
if (nskip > 0) { nskip. <- 1:nskip } else { nskip. <- 0 }
if (keepevery > 1) {
  res_ctrl$sigma_rwd <- c(res_ctrl$sigma_rwd[nskip.],
                           res_ctrl$sigma_rwd[nskip + seq(1, ndpost * keepevery, keepevery)])
}
res_ctrl$sigma_rwd <- res_ctrl$sigma_rwd[-(nskip.)]

cat(sprintf("  Control model done. sigma (shared): mean=%.3f\n",
            mean(res_ctrl$sigma_rwd)))

# ============================================================
# 6. TREATMENT ARM: standard BART (aBART) on EloKRD
# ============================================================
cat("\n=== Fitting TREATMENT arm model (aBART on EloKRD) ===\n")

x.train.trt <- data.frame(X = as.matrix(X_all[trt_idx, ]))
y.train.trt <- log(Y[trt_idx])
event.trt   <- as.integer(event[trt_idx])
x.test.trt  <- data.frame(X = as.matrix(X_all[trt_idx, ]))  # predict on same patients

xinfo2 <- matrix(0, 0, 0)
temp2 <- bartModelMatrix(x.train.trt, numcut, usequants = FALSE,
                         xinfo = xinfo2, rm.const = TRUE)
x.train.trt.m <- t(temp2$X)
numcut.trt    <- temp2$numcut
xinfo.trt     <- temp2$xinfo

x.test.trt.m <- bartModelMatrix(x.test.trt)
x.test.trt.m <- t(x.test.trt.m[, temp2$rm.const])

p_trt  <- nrow(x.train.trt.m)
n_trt  <- ncol(x.train.trt.m)
np_trt <- ncol(x.test.trt.m)

cat(sprintf("  Training: n=%d, p=%d; Test: np=%d\n", n_trt, p_trt, np_trt))

# Center on observed mean so BART learns residuals, not absolute level
trt_offset <- mean(y.train.trt)
# trt_offset <- 0   # uncentered version (commented out)
y.trt.scaled <- y.train.trt - trt_offset
cat(sprintf("  Centering offset: %.3f (exp=%.2f yr)\n", trt_offset, exp(trt_offset)))

df_trt <- data.frame(t(x.train.trt.m))
aft_fit_trt <- survreg(Surv(exp(y.trt.scaled), event.trt) ~ .,
                       data = df_trt, dist = "lognormal")
sigest_trt <- aft_fit_trt$scale

qchi <- qchisq(1.0 - sigquant, sigdf)
lambda_trt <- (sigest_trt^2 * qchi) / sigdf
# Treatment arm uses original BART hyperparameters (ntree=50, k=2);
# control arm is the one that was shrunk (ntree=5, k=10) to suppress
# per-draw tree-prior swings in the EloKRd counterfactual prediction.
ntree_trt <- 50L
k_trt     <- 2
tau_trt   <- (max(y.trt.scaled) - min(y.trt.scaled)) / (2 * k_trt * sqrt(ntree_trt))

ntype <- 1
sparse <- FALSE
theta <- 0
omega <- 1
a <- 0.5
b <- 1
augment <- FALSE
rho <- p_trt
w <- rep(1, n_trt)
grp <- temp2$grp
if (is.null(grp)) grp <- 1:p_trt

cat("  Running aBART MCMC...\n")
res_trt <- cabart(ntype,
                  n_trt,
                  p_trt,
                  np_trt,
                  x.train.trt.m,
                  y.trt.scaled,
                  event.trt,
                  x.test.trt.m,
                  ntree_trt,
                  numcut.trt,
                  ndpost * keepevery,
                  nskip,
                  keepevery,
                  beta_tree,
                  alpha,
                  trt_offset,  # offset (centering)
                  tau_trt,
                  sigdf,
                  lambda_trt,
                  sigest_trt,
                  w,
                  sparse,
                  theta,
                  omega,
                  grp,
                  a,
                  b,
                  rho,
                  augment,
                  10000L,  # printevery
                  xinfo.trt,
                  42L)

if (nskip > 0) { nskip. <- 1:nskip } else { nskip. <- 0 }
if (keepevery > 1) {
  res_trt$sigma <- c(res_trt$sigma[nskip.],
                     res_trt$sigma[nskip + seq(1, ndpost * keepevery, keepevery)])
}
res_trt$sigma <- res_trt$sigma[-(nskip.)]

cat(sprintf("  Treatment model done. sigma: mean=%.3f\n", mean(res_trt$sigma)))

# ============================================================
# 7. COMPUTE POPULATION RMST (restricted mean survival time)
# ============================================================
# Primary estimand: RMST ratio at tau = 5 years.  tau is chosen
# inside the observed follow-up horizon so the integral is
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

rmst_ctrl <- compute_pop_rmst(res_ctrl$f.test,  res_ctrl$sigma_rwd, tau_rmst,
                              "Control (UCMM->EloKRD via theta_1, shared sigma)")
rmst_trt  <- compute_pop_rmst(res_trt$yhat.test, res_trt$sigma,     tau_rmst,
                              "Treatment (EloKRD)")

# ============================================================
# 8. POSTERIOR TREATMENT EFFECT
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
# 9b. ARM-SPECIFIC ESTIMATES (population RMST + residual sigma, median [95% CI])
# ============================================================
arm_summ <- function(x) {
  q <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
  c(est = median(x, na.rm = TRUE), lo = unname(q[1]), hi = unname(q[2]))
}
rmst_trt_est   <- arm_summ(rmst_trt)
rmst_ctrl_est  <- arm_summ(rmst_ctrl)
sigma_trt_est  <- arm_summ(res_trt$sigma)
sigma_ctrl_est <- arm_summ(res_ctrl$sigma_rwd)
cat("\n=== Arm-specific estimates (median [95% CI]) ===\n")
cat(sprintf("  RMST  treatment: %.3f  [%.3f, %.3f]\n", rmst_trt_est["est"],  rmst_trt_est["lo"],  rmst_trt_est["hi"]))
cat(sprintf("  RMST  control  : %.3f  [%.3f, %.3f]\n", rmst_ctrl_est["est"], rmst_ctrl_est["lo"], rmst_ctrl_est["hi"]))
cat(sprintf("  sigma treatment: %.3f  [%.3f, %.3f]\n", sigma_trt_est["est"],  sigma_trt_est["lo"],  sigma_trt_est["hi"]))
cat(sprintf("  sigma control  : %.3f  [%.3f, %.3f]\n", sigma_ctrl_est["est"], sigma_ctrl_est["lo"], sigma_ctrl_est["hi"]))

# ============================================================
# 10. SAVE RESULTS
# ============================================================
out_dir <- file.path(projDir, "res")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results <- list(
  post_ratio   = post_ratio_valid,
  delta_hat    = delta_hat,
  ci_95        = ci_95,
  delta_diff   = delta_diff,
  ci_diff_95   = ci_diff_95,
  rmst_ctrl    = rmst_ctrl,
  rmst_trt     = rmst_trt,
  rmst_trt_est   = rmst_trt_est,
  rmst_ctrl_est  = rmst_ctrl_est,
  sigma_trt_est  = sigma_trt_est,
  sigma_ctrl_est = sigma_ctrl_est,
  tau_rmst     = tau_rmst,
  res_ctrl     = res_ctrl,
  res_trt      = res_trt,
  settings     = list(ndpost = ndpost, nskip = nskip, ntree = ntree, k = k,
                      prior = prior_val, alpha = alpha, beta = beta_tree,
                      target_N = target_N)
)

out_file <- file.path(out_dir,
                       sprintf("mBART_results_%s_%s_nt%d_k%g_N%d.RData",
                               OUTCOME, data_tag, ntree, k, target_N))
saveRDS(results, file = out_file)
cat(sprintf("\nResults saved to: %s\n", out_file))

mbart_summary[[length(mbart_summary) + 1]] <- data.frame(
  ntree = ntree, k = k, target_N = target_N,
  delta_hat = delta_hat, ci_lo = unname(ci_95[1]), ci_hi = unname(ci_95[2]),
  row.names = NULL)

}  # end loop over N_target multipliers
}  # end loop over control-arm configs

# =====================================================================
# Control-arm sensitivity comparison: RMST(tau) ratio (EloKRd/UCMM)
# per (ntree, k) and N_target.
# =====================================================================
cat("\n", strrep("=", 64),
    "\n=== mBART control-arm sensitivity (RMST ratio EloKRd/UCMM) ===\n",
    strrep("=", 64), "\n", sep = "")
if (length(mbart_summary)) {
  mbart_tab <- do.call(rbind, mbart_summary)
  print(transform(mbart_tab,
                  delta_hat = round(delta_hat, 3),
                  ci_lo     = round(ci_lo, 3),
                  ci_hi     = round(ci_hi, 3)),
        row.names = FALSE)
} else {
  cat("No MAP-BART runs (no config had a matching ess calibration).\n")
}

cat("\n=== DONE ===\n")
