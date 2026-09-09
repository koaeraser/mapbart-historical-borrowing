## ess_calibration.R (mapbart-case-study-mm copy)
##
## Prior ESS calibration for MAP-AFT-BART on the mapbart-case-study-mm merged EloKRD /
## UCMM dataset.  Duplicate of mapbart-sim-survival/ess_local/ess_calibration.R
## with two changes:
##
##   * Canonical RWD: data_cleaned/merged_elokrd_ucmm_n283.RData (default),
##     restricted to source == "UCMM" (the external-control arm; 253 of the
##     283 subjects in the n283 file, 48 of 78 in the n78 file).
##     Survival outcome PFS (default) or OS via --outcome.
##   * Covariates: age, sex, race, ethnicity, high_risk_cyto,
##     transplant_off_protocol (same set as mapbart-case-study-mm/mBART_analysis.R;
##     categorical columns are one-hot encoded).
##
## Infrastructure (cabart_ess.cpp, include_ess/, ess_plot.R) is shared
## via symlinks pointing at mapbart-sim-survival/ess_local/.  RBesT_lite/ is
## also a symlink, pointing at mapbart-sim-gaussian/ess_local/RBesT_lite/ (the
## canonical Stan/gMAP copy).  Outputs go to ./res/ under this folder.
##
## Usage:
##   Rscript ess_calibration.R [--data <RData file>]
##                                   [--source UCMM]
##                                   [--outcome OS|PFS]
##                                   [--s2-from 0.1] [--s2-to 0.9] [--s2-by 0.025]
##                                   [--s2-grid "..."]   # explicit s^2 values; overrides from/to/by (empty default -> use from/to/by)
##                                   [--Ntarget 30]  [--seed 6]
##                                   [--map-model synthetic] [--beta-prior auto] [--sigma-ref auto]
##                                   [--ntree 50] [--k 2]   # AFT-BART control config; tags output _nt<ntree>_k<k>
##                                   [--B 100] [--burn 1000] [--q 0.95]
##
## See the parent script in mapbart-sim-survival/ess_local/ for full flag
## documentation; this copy keeps the same CLI surface.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppEigen)
  library(parallel)
  library(survival)
})

scriptDir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)),
  error = function(e) "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/ess_local"
)
if (!nzchar(scriptDir) || is.na(scriptDir))
  scriptDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/ess_local"

## Shared ess_local infrastructure now lives in this folder via
## symlinks (cabart_ess.cpp, include_ess/, RBesT_lite/, ess_plot.R),
## so it's all addressable through scriptDir like the other ess_local
## folders.  bartModelMatrix.R is the one absolute reference left.
CABART_ESS   <- file.path(scriptDir, "cabart_ess.cpp")
RBESTLITE    <- file.path(scriptDir, "RBesT_lite", "load.R")
BARTMODELMAT <- "/Users/oliviazhang/Desktop/bartModelMatrix.R"

## --- CLI ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  k <- grep(paste0("^", flag, "$"), args)
  if (length(k) == 0) return(default)
  args[k + 1]
}

# Pass --data as just the FILE NAME; it is resolved against the project's
# data_cleaned/ directory (a full path also works). --source defaults to UCMM.
# The cohort is identified by the data file name -- its n<N> token tags the
# output, which mBART_analysis.R matches by file name.
clean_dir     <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/data_cleaned"
default_data  <- file.path(clean_dir, "merged_elokrd_ucmm_n283.RData")
data_path     <- get_arg("--data",        default_data)
# Resolve a bare filename (or any not-found path) against data_cleaned/.
if (!file.exists(data_path) && file.exists(file.path(clean_dir, basename(data_path))))
  data_path <- file.path(clean_dir, basename(data_path))
source_filter <- get_arg("--source",      "UCMM")
outcome       <- toupper(get_arg("--outcome", "PFS"))
stopifnot(outcome %in% c("OS", "PFS"))
s2_from       <- as.numeric(get_arg("--s2-from",   0.1))
s2_to         <- as.numeric(get_arg("--s2-to",     0.9))
s2_by         <- as.numeric(get_arg("--s2-by",     0.025))
N_target      <- as.numeric(get_arg("--Ntarget",   30))
seed          <- as.integer(get_arg("--seed",      6L))
map_model  <- get_arg("--map-model", "synthetic")
## --beta-prior is the SD of the theta_h prior in standard/synthetic.  Default
## "auto" sets it to the BART leaf-mean prior SD tau_theta (= itau),
## so the MAP anchor prior N(0, beta.prior^2) on the leaf mean matches
## the prior BART actually used to generate that leaf value.  Resolved
## after tau_theta is computed in Stage 1.
beta_prior_arg <- get_arg("--beta-prior", "auto")
sigma_ref_arg <- get_arg("--sigma-ref", "auto")
## --ntree / --k: the AFT-BART control config the calibration is run under
## (should match the MAP-BART control config it will calibrate, since tau_theta =
## (max-min)/(2*k*sqrt(ntree)) and the leaf partial-residual structure depend on
## both).  The config tag _nt<ntree>_k<k> is appended to every output filename so
## calibrations for different configs never clobber each other; mBART_analysis.R
## matches the calibration by this tag.
ntree_arg <- as.integer(get_arg("--ntree", 50L))
k_arg     <- as.numeric(get_arg("--k",     2))
cfg_tag   <- sprintf("_nt%d_k%g", ntree_arg, k_arg)
cat("BART config tag:", sub("^_", "", cfg_tag), "\n")
B             <- as.integer(get_arg("--B",         100L))
burn          <- as.integer(get_arg("--burn",      1000L))
q_target      <- as.numeric(get_arg("--q",         0.95))
stopifnot(map_model %in% c("standard", "hybrid", "synthetic"))

tilde_nu <- 3
## --s2-grid: explicit comma/space-separated s^2 values to sweep.  Empty by
## default -- in which case the grid is built from --s2-from/--s2-to/--s2-by
## (default 0.1..0.9 by 0.025, 33 pts), so there is no opaque baked-in grid.
## Pass --s2-grid "v1,v2,..." to override with a specific (e.g. re-tuned) set.
s2_grid_arg <- get_arg("--s2-grid", "")
if (nzchar(s2_grid_arg)) {
  S_grid <- sort(unique(as.numeric(strsplit(s2_grid_arg, "[, ]+")[[1]])))
} else {
  S_grid <- seq(s2_from, s2_to, by = s2_by)
}

cat(sprintf("ess_calibration (mapbart-case-study-mm): source=%s  outcome=%s  data=%s  B=%d  burn=%d  N_target=%g  q=%g  gMAP=%s\n",
            source_filter, outcome, basename(data_path), B, burn, N_target, q_target, map_model))
cat(sprintf("S_grid = [%s]  (K=%d points)\n",
            paste(S_grid, collapse=", "), length(S_grid)))

## --- canonical RWD ------------------------------------------------------
## File contains a single data.frame `merged`; load via load() not readRDS.

if (!file.exists(data_path))
  stop("--data path not found: '", data_path,
       "'  (omit --data to use the default: ", default_data, ")")
cat("Canonical RWD file:", basename(data_path), "\n")
.env <- new.env()
load(data_path, envir = .env)
stopifnot("merged" %in% ls(.env))
merged <- .env$merged
stopifnot("source" %in% colnames(merged))

## Output tag = the n<N> token from the data file name (falls back to the
## file stem). mBART_analysis.R matches the calibration file by this tag.
cohort_tag <- regmatches(basename(data_path), regexpr("n[0-9]+", basename(data_path)))
if (!length(cohort_tag)) cohort_tag <- sub("\\.RData$", "", basename(data_path))
cat(sprintf("Data file tag: %s\n", cohort_tag))

## Reproducibility hygiene: clear THIS run's stale output(s) up front, so a
## crash can never leave a previous result on disk to be mistaken for fresh
## output.  Scoped to the file(s) this invocation writes (this source +
## outcome + map_model + cohort_tag) -- sibling calibrations are untouched.
.resDir <- file.path(scriptDir, "res")
.stale  <- file.path(.resDir, c(
  sprintf("ess_calibration_%s_%s_%s_%s%s.RData",
          tolower(source_filter), tolower(outcome), map_model, cohort_tag, cfg_tag),
  sprintf("ess_calibration_%s_%s_%s_%s%s_checkpoint.RData",
          tolower(source_filter), tolower(outcome), map_model, cohort_tag, cfg_tag)))
.stale  <- .stale[file.exists(.stale)]
if (length(.stale)) {
  cat("Clearing stale output before rerun:\n    ",
      paste(basename(.stale), collapse = "\n     "), "\n", sep = "")
  invisible(file.remove(.stale))
}

merged <- merged[merged$source == source_filter, , drop = FALSE]
cat(sprintf("After filter source==\"%s\": n=%d\n",
            source_filter, nrow(merged)))

## Outcome columns.
y_col <- if (outcome == "OS") "os_months"  else "pfs_months"
e_col <- if (outcome == "OS") "os_status"  else "pfs_status"
stopifnot(y_col %in% colnames(merged), e_col %in% colnames(merged))
Yobs_RWD  <- as.numeric(merged[[y_col]])
delta_RWD <- as.integer(merged[[e_col]])

## Covariate set: matches mapbart-case-study-mm/mBART_analysis.R.
covariate_cols <- c("age", "sex", "race", "ethnicity",
                     "high_risk_cyto", "transplant_off_protocol")
prep_covariates <- function(df, cols) {
  X <- data.frame(row.names = seq_len(nrow(df)))
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
        for (k in 2:length(lvls))
          X[[paste0(v, "_", gsub("[/ ]", "_", lvls[k]))]] <-
            as.integer(vals == lvls[k])
      }
    }
  }
  X
}
X_RWD <- as.matrix(prep_covariates(merged, covariate_cols))
n_RWD <- length(Yobs_RWD)
p     <- ncol(X_RWD)
log_Y <- log(pmax(Yobs_RWD, .Machine$double.xmin))
cat(sprintf("RWD subset: n=%d, p=%d, events=%d (%.0f%%)\n",
            n_RWD, p, sum(delta_RWD), 100 * mean(delta_RWD)))
cat("Covariate columns:", paste(colnames(X_RWD), collapse = ", "), "\n")

## --- Stage 1: leaf summaries -------------------------------------------

H        <- ntree_arg     # AFT-BART ntree (--ntree; default 50)
alpha    <- 0.95
beta     <- 2
k_tau    <- k_arg         # leaf-prior scale (--k; default 2): tau_theta = (max-min)/(2*k*sqrt(H))
numcut   <- 100L
nu       <- 3
sigquant <- 0.90

.aft_fit  <- survreg(Surv(Yobs_RWD, delta_RWD) ~ ., data = data.frame(X_RWD),
                      dist = "lognormal")
sigest    <- .aft_fit$scale
lambda    <- (sigest^2 * qchisq(1 - sigquant, nu)) / nu
tau_theta <- (max(log_Y) - min(log_Y)) / (2 * k_tau * sqrt(H))
cat(sprintf("sigest(AFT scale)=%.3f  lambda=%.3f  tau_theta=%.3f\n",
            sigest, lambda, tau_theta))

## Resolve --beta-prior now that tau_theta (the BART leaf-mean prior SD,
## passed to cabart_ess as itau) is known.  "auto" => beta.prior =
## tau_theta so the standard/synthetic theta_h prior N(0, beta.prior^2) matches
## the BART leaf prior; a numeric overrides.
beta_prior <- if (identical(beta_prior_arg, "auto")) tau_theta else
                as.numeric(beta_prior_arg)
cat(sprintf("beta.prior (theta_h prior SD) = %.4f  (%s)\n", beta_prior,
            if (identical(beta_prior_arg, "auto"))
              "auto: BART leaf-prior SD tau_theta" else "user-specified"))

source(BARTMODELMAT)
temp    <- bartModelMatrix(X_RWD, numcut, usequants = FALSE, cont = FALSE,
                            xinfo = matrix(0, 0, 0), rm.const = TRUE)
Xt      <- t(temp$X)
numcut  <- temp$numcut
xinfo   <- temp$xinfo
p_eff   <- nrow(Xt)
cat(sprintf("After bartModelMatrix: p_eff = %d (rm.const removed %d)\n",
            p_eff, p - p_eff))

t0 <- Sys.time()
cat("Compiling", CABART_ESS, "... ")
sourceCpp(CABART_ESS, verbose = FALSE)
cat(sprintf("done in %s\n", format(Sys.time() - t0)))

cat(sprintf("Running AFT-BART on RWD (H=%d, B=%d, burn=%d) ...\n",
            H, B, burn))
set.seed(seed)
t0 <- Sys.time()
bart_res <- cabart_ess(
  in_       = ncol(Xt),
  ip        = nrow(Xt),
  ix        = as.numeric(Xt),
  iy        = as.numeric(log_Y),
  idelta    = as.integer(delta_RWD),
  im        = as.integer(H),
  inc       = as.integer(numcut),
  ind       = as.integer(B),
  iburn     = as.integer(burn),
  ipower    = as.numeric(beta),
  ibase     = as.numeric(alpha),
  itau      = as.numeric(tau_theta),
  inu       = as.numeric(nu),
  ilambda   = as.numeric(lambda),
  isigest   = as.numeric(sigest),
  iiw       = rep(1.0, ncol(Xt)),
  iverbose  = 1L,
  iXinfo    = if (is.matrix(xinfo) && nrow(xinfo) > 0) xinfo else NULL
)
cat(sprintf("Stage 1 done in %s.  leaf_table: %d rows\n",
            format(Sys.time() - t0), nrow(bart_res$leaf_table)))

leaf_table <- as.data.frame(bart_res$leaf_table)

use_auto_sigma <- identical(sigma_ref_arg, "auto")
sigma_ref      <- if (use_auto_sigma) NA_real_ else as.numeric(sigma_ref_arg)
cat(sprintf("sigma_ref = %s\n",
            if (use_auto_sigma) "auto: per-draw BART sigma_b (AFT log-scale)"
            else sprintf("%.4f (user-specified)", sigma_ref)))

## --- Stage 2: RBesT_lite gMAP + ELIR ESS -------------------------------

source(RBESTLITE)

options(
  RBesT.MC.chains  = 1,
  RBesT.MC.warmup  = 500,
  RBesT.MC.iter    = 2000,
  RBesT.MC.thin    = 1,
  RBesT.MC.control = list(adapt_delta = 0.95,
                          max_treedepth = 10,
                          stepsize = 0.1),
  mc.cores         = 1L
)
AUTOMIX_NC     <- 1:3
AUTOMIX_THRESH <- 0

ess_one_row <- function(ybar, n_leaf, sigma_b, S_grid, tilde_nu,
                          map_model, beta_prior, sigma_for_ess) {
  SE <- sigma_b / sqrt(n_leaf)
  one_row <- data.frame(y.bar = ybar, n = n_leaf, y.se = SE,
                         study = factor(1, levels = "1"))
  ess_vec <- numeric(length(S_grid))
  for (k in seq_along(S_grid)) {
    s2 <- S_grid[k]
    args_gmap <- list(
      formula     = cbind(y.bar, y.se) ~ 1 | study,
      data        = one_row,
      tau.dist    = "InvGamma",
      tau.prior   = c(tilde_nu / 2, tilde_nu * s2 / 2),
      family      = gaussian,
      version     = map_model
    )
    if (map_model %in% c("standard", "synthetic"))
      args_gmap$beta.prior <- beta_prior
    MAP_prior   <- do.call(gMAP, args_gmap)
    MAP_mixture <- automixfit(MAP_prior, Nc = AUTOMIX_NC,
                               thresh = AUTOMIX_THRESH)
    ess_vec[k]  <- ess(MAP_mixture, sigma = sigma_for_ess)
  }
  ess_vec
}

process_b <- function(b, leaf_table, S_grid, tilde_nu, H, N_target,
                       map_model, beta_prior, sigma_ref, seed_base) {
  ## Deterministic per-draw seed: makes this b's gMAP fits (whose Stan
  ## seeds are drawn from R's RNG) depend only on (seed_base, b), not on
  ## which worker runs it, the load-balancing order, or the core count.
  set.seed(seed_base + b)
  t_start <- Sys.time()
  rows_b   <- leaf_table[leaf_table$b == b, ]
  K        <- length(S_grid)
  ESS_h    <- matrix(0, nrow = H, ncol = K)
  n_rows   <- nrow(rows_b)
  message(sprintf("[%s] b=%-3d start  (%d leaves x %d s^2 = %d gMAP fits, pid=%d)",
                  format(t_start, "%H:%M:%S"), b, n_rows, K,
                  n_rows * K, Sys.getpid()))
  for (r in seq_len(n_rows)) {
    sigma_for_ess <- if (is.na(sigma_ref)) rows_b$sigma[r] else sigma_ref
    ess_vec <- tryCatch(
      ess_one_row(rows_b$ybar[r], rows_b$n[r], rows_b$sigma[r],
                   S_grid, tilde_nu, map_model, beta_prior,
                   sigma_for_ess),
      error = function(e) {
        message(sprintf("    b=%d h=%d leaf=%d FAILED: %s",
                        b, rows_b$h[r], rows_b$leaf[r],
                        conditionMessage(e)))
        rep(NA_real_, K)
      })
    h_idx <- rows_b$h[r] + 1
    ESS_h[h_idx, ] <- ESS_h[h_idx, ] + ess_vec
  }
  ESSbar_b <- apply(ESS_h, 2, median, na.rm = TRUE)  # median over H trees
                                                     # (robust to tree-level
                                                     #  ESS spikes; cf. mean
                                                     #  in Algorithm 1)
  dt <- as.numeric(Sys.time() - t_start, units = "secs")
  message(sprintf("[%s] b=%-3d done   (%.1fs)  ESSbar range [%.1f, %.1f]",
                  format(Sys.time(), "%H:%M:%S"), b, dt,
                  min(ESSbar_b, na.rm = TRUE),
                  max(ESSbar_b, na.rm = TRUE)))
  list(b = b, ESS_h = ESS_h, ESSbar_b = ESSbar_b,
        below = as.integer(ESSbar_b <= N_target))
}

cores <- as.integer(Sys.getenv("CORES",
                                unset = max(1L, parallel::detectCores() - 1L)))
cells_per_b <- sum(leaf_table$b == 0) * length(S_grid)
total_cells <- nrow(leaf_table) * length(S_grid)
cat(sprintf("\nStage 2: %d cores, ~%d Stan fits per b, B=%d work items, %d total fits\n",
            cores, cells_per_b, B, total_cells))
flush.console()

## PSOCK cluster, NOT mclapply: gMAP calls rstan::sampling, and
## rstan/RcppParallel/TBB thread state does not survive a fork() -- under
## mclapply the Stan sampler intermittently fails inside forked children,
## silently NA-ing a large fraction of draws.  PSOCK workers are fresh R
## processes; outfile="" forwards their progress to this terminal.
## CORES=1 runs sequentially in-process.  RBesT_lite is loaded from the
## shared ess_local copy in each worker.
b_seq <- seq_len(B) - 1L
rbeslite_dir <- dirname(RBESTLITE)
gmap_opts <- list(
  RBesT.MC.chains  = 1, RBesT.MC.warmup = 500, RBesT.MC.iter = 2000,
  RBesT.MC.thin    = 1,
  RBesT.MC.control = list(adapt_delta = 0.95, max_treedepth = 10,
                          stepsize = 0.1),
  mc.cores         = 1L)

t0 <- Sys.time()
if (cores <= 1L) {
  worker_results <- lapply(
    b_seq, process_b,
    leaf_table = leaf_table, S_grid = S_grid, tilde_nu = tilde_nu, H = H,
    N_target = N_target, map_model = map_model,
    beta_prior = beta_prior, sigma_ref = sigma_ref, seed_base = seed)
} else {
  cl <- parallel::makeCluster(cores, type = "PSOCK", outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl,
    varlist = c("ess_one_row", "process_b", "rbeslite_dir", "gmap_opts"),
    envir = environment())
  parallel::clusterEvalQ(cl, {
    Sys.setenv(RBESLITE_DIR = rbeslite_dir)
    source(file.path(rbeslite_dir, "load.R"))
    options(gmap_opts)
    AUTOMIX_NC     <- 1:3
    AUTOMIX_THRESH <- 0
    NULL
  })
  ## No clusterSetRNGStream needed: process_b seeds itself per-b from
  ## seed_base, so results are independent of worker assignment and the
  ## load-balanced dispatch order.
  worker_results <- parallel::parLapplyLB(
    cl, b_seq, process_b,
    leaf_table = leaf_table, S_grid = S_grid, tilde_nu = tilde_nu, H = H,
    N_target = N_target, map_model = map_model,
    beta_prior = beta_prior, sigma_ref = sigma_ref, seed_base = seed)
}
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\nStage 2 done in %.1fs (%.1fm)\n", elapsed, elapsed / 60))

## Pre-aggregate checkpoint: dump raw worker_results to disk *before*
## the aggregation block.  If aggregation crashes (e.g. a dead worker
## returned a non-list), the raw output is preserved and aggregation
## can be retried from disk without re-running Stage 2.
resDir <- file.path(scriptDir, "res")
dir.create(resDir, showWarnings = FALSE)
.ckpt_path <- file.path(resDir,
        sprintf("ess_calibration_%s_%s_%s_%s%s_checkpoint.RData",
                tolower(source_filter), tolower(outcome), map_model, cohort_tag, cfg_tag))
saveRDS(list(worker_results = worker_results,
              leaf_table     = leaf_table,
              S_grid         = S_grid,
              N_target       = N_target,
              q              = q_target,
              B              = B,
              sigma_draws    = bart_res$sigma_draws),
         .ckpt_path)
cat("Stage-2 checkpoint:", .ckpt_path, "\n")

## --- Aggregate + selection rule ----------------------------------------
## Defensive: a worker may have died mid-run (Stan crash, OOM, etc.) and
## returned a try-error or NULL instead of a list.  Drop those, warn,
## and aggregate over the survivors.

K          <- length(S_grid)
ok_workers <- vapply(worker_results, function(w) {
  is.list(w) && !is.null(w$ESSbar_b) && is.numeric(w$ESSbar_b) &&
    length(w$ESSbar_b) == K
}, logical(1))
n_bad <- sum(!ok_workers)
if (n_bad > 0) {
  bad_idx <- which(!ok_workers)
  cat(sprintf("WARNING: %d/%d workers returned no valid ESSbar_b (b = %s).\n",
              n_bad, length(worker_results),
              paste(bad_idx - 1L, collapse = ", ")))
  cat("  Aggregation will proceed over the surviving workers.\n")
}
valid_results <- worker_results[ok_workers]
B_used <- length(valid_results)
stopifnot(B_used > 0)

ESSbar <- t(sapply(valid_results, `[[`, "ESSbar_b"))   # B_used x K

## Derive c_s2 / Ghat directly from ESSbar with na.rm so a draw whose
## ESSbar_b is NA in one s^2 column (every tree had a failed-gMAP leaf
## there) does not poison the whole column.  Each column's Ghat is the
## fraction of *valid* draws with ESSbar_b <= N_target.
below_mat <- ESSbar <= N_target                        # may contain NA
c_s2      <- colSums(below_mat, na.rm = TRUE)
n_valid   <- colSums(!is.na(ESSbar))                   # valid draws per s^2
Ghat      <- ifelse(n_valid > 0, c_s2 / n_valid, NA_real_)

satisfies   <- which(Ghat >= q_target)
hat_s2_star <- if (length(satisfies)) S_grid[min(satisfies)] else NA_real_

cat("\n--- Calibration results ---\n")
print(data.frame(
  s2          = S_grid,
  Ghat        = round(Ghat, 3),
  c_s2        = c_s2,
  n_valid     = n_valid,
  medH_meanB_ESS = round(colMeans(ESSbar, na.rm = TRUE), 1)
))
if (is.na(hat_s2_star)) {
  cat(sprintf("\nNo s^2 in grid satisfies Pr{ESSbar <= %g} >= %g.\n",
              N_target, q_target),
       "  Extend grid upward and rerun.\n")
} else {
  cat(sprintf("\nCalibrated s^2_star = %.3f  (smallest with Ghat >= %g)\n",
              hat_s2_star, q_target))
}

## --- Save ---------------------------------------------------------------

out <- list(
  canonical_file = data_path,
  source_filter  = source_filter,
  outcome        = outcome,
  n_RWD          = n_RWD,
  p              = p_eff,
  covariates     = colnames(X_RWD),
  S_grid         = S_grid,
  Ghat           = Ghat,
  c_s2           = c_s2,
  ESSbar         = ESSbar,
  sigma_draws    = bart_res$sigma_draws,
  leaf_table     = leaf_table,
  N_target       = N_target,
  q              = q_target,
  hat_s2_star    = hat_s2_star,
  B              = B,
  B_used         = B_used,
  burn           = burn,
  H              = H,
  ntree          = ntree_arg,
  k              = k_arg,
  tilde_nu       = tilde_nu,
  map_model   = map_model,
  beta_prior     = beta_prior,
  sigma_ref      = sigma_ref
)
out_file <- file.path(resDir,
        sprintf("ess_calibration_%s_%s_%s_%s%s.RData",
                tolower(source_filter), tolower(outcome), map_model, cohort_tag, cfg_tag))
saveRDS(out, file = out_file)
cat("Wrote", out_file, "\n")
