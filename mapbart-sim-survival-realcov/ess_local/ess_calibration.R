## ess_calibration.R -- Prior ESS calibration for MAP-AFT-BART
## (censored log-normal outcome) via AFT-BART partial residuals.
## Survival counterpart of ../../mapbart-sim-gaussian-realcov/ess_local/ess_calibration.R.
##
## There is no vanilla-MAP survival companion (no MAP_surv.R analysis
## script exists in mapbart-sim-survival/), so this script runs alone --
## ess_plot.R --outcome surv shows AFT-BART boxes without a vanilla-MAP
## overlay.
##
## Stage 1 (C++): cabart_ess.cpp runs AFT-BART (heterbart + per-sweep
##   rtnorm imputation of censored latent log-times) on the canonical
##   RWD and emits a leaf_table with rows (b, h, leaf, ybar, n, sigma).
##   ybar is the within-leaf mean of the partial residual on the latent
##   log-time scale at the moment tree h is being updated; sigma is the
##   AFT log-scale SD in scope at that point (sigma^(b-1)).
##
## Stage 2 (R): identical to the Gaussian script -- for each leaf_table
##   row x each s^2 in the candidate grid, fit the MAP model with
##   RBesT_lite::gMAP (family = gaussian on the log-time scale),
##   approximate the predictive by automixfit, and compute ELIR ESS via
##   ess().  Per draw b: sum ESS across the leaves of each tree h, then
##   take the median over the H trees -> ESSbar_b (robust to tree-level
##   ESS spikes).  Apply the Pr-rule
##     Ghat(s^2) = #{ b : ESSbar_b <= N_target } / B
##     hat_s2_star = min{ s^2 in grid : Ghat(s^2) >= q }.
##
## Usage:
##   Rscript ess_calibration.R [--data <RWD file or dir>]
##                                   [--s2-from 0.15] [--s2-to 0.55] [--s2-by 0.05]
##                                   [--s2-grid "..."]   # per-target clusters (default); overrides from/to/by
##                                   [--Ntarget 100] [--seed 6]
##                                   [--map-model hybrid] [--outcome n200]
##                                   [--beta-prior auto] [--sigma-ref auto]
##                                   [--B 100] [--burn 1000] [--q 0.95]
##
## CLI is identical to the Gaussian ess_calibration.R; only the data
## and Stage-1 machinery (cabart_ess vs cwbart_ess) differ.  --sigma-ref controls sigma passed to ess(...):
##   auto (default)  -- per-draw BART sigma_b (AFT log-scale SD)
##   <numeric>       -- fixed value for every fit.
##
## Defaults: canonical RWD = first *.RData file (alphabetical) under
## /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-realcov/data_v3/, D == 0 & Z == 0
## subset.  Also reads env var CORES (defaults to detectCores() - 1)
## for the per-b parallel fan-out in Stage 2.

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppEigen)
  library(parallel)
  library(survival)
})

scriptDir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)),
  error = function(e) "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-realcov/ess_local"
)
if (!nzchar(scriptDir) || is.na(scriptDir))
  scriptDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-realcov/ess_local"

## --- CLI ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  k <- grep(paste0("^", flag, "$"), args)
  if (length(k) == 0) return(default)
  args[k + 1]
}

default_data  <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-realcov/data_v3"
data_path     <- get_arg("--data",        default_data)
s2_from       <- as.numeric(get_arg("--s2-from",   0.15))
s2_to         <- as.numeric(get_arg("--s2-to",     0.55))
s2_by         <- as.numeric(get_arg("--s2-by",     0.05))
N_target      <- as.numeric(get_arg("--Ntarget",   100))
seed          <- as.integer(get_arg("--seed",      6L))
map_model  <- get_arg("--map-model", "synthetic")
## --outcome: optional RWD-size label (e.g. n200 / n30) appended to the
## output filename so calibrations for different RWD sizes don't overwrite
## each other; pass the matching file to ess_plot.R via --bart.
rwd_outcome   <- get_arg("--outcome", "")
out_tag       <- if (nzchar(rwd_outcome)) paste0("_", rwd_outcome) else ""
## Scenario tag for res/ filenames (sc1E / sc2E / sc3E_cor<rho>), parsed from the
## source RWD file so different scenarios never clobber each other.
.sc_src <- if (dir.exists(data_path)) {
  .fs <- sort(list.files(data_path, pattern = "RData$")); if (length(.fs)) .fs[1] else ""
} else basename(data_path)
.m        <- regmatches(.sc_src, regexpr("sc[0-9]+[A-Za-z]*(_cor-?[0-9.]+)?", .sc_src))
sc_tag    <- if (length(.m)) .m else ""
sc_suffix <- if (nzchar(sc_tag)) paste0("_", sc_tag) else ""
if (nzchar(sc_tag)) cat("Scenario tag:", sc_tag, "\n")
## --beta-prior is the SD of the theta_h/mu prior in standard/synthetic (ignored
## by hybrid).  "auto" sets it to the BART leaf-mean prior SD
## tau_theta (= itau), so the MAP prior on the leaf location matches the
## prior AFT-BART used to generate it.  Resolved after tau_theta below.
beta_prior_arg <- get_arg("--beta-prior", "auto")
sigma_ref_arg <- get_arg("--sigma-ref", "auto")
## --ntree / --k: the CONTROL-arm AFT-BART config the calibration is run under
## (should match the MAP-BART control config it will calibrate, since tau_theta =
## (max-min)/(2*k*sqrt(ntree)) and the leaf partial-residual structure depend on
## both).  Defaults (50, 2) reproduce the canonical config and the UNtagged
## filename; any other config tags the output _nt<ntree>_k<k>.
ntree_arg <- as.integer(get_arg("--ntree", 50L))
k_arg     <- as.numeric(get_arg("--k",     2))
cfg_tag   <- sprintf("_nt%d_k%g", ntree_arg, k_arg)
cat("BART config tag:", sub("^_", "", cfg_tag), "\n")
## --- BART-only flags ---
B             <- as.integer(get_arg("--B",         100L))
burn          <- as.integer(get_arg("--burn",      1000L))
q_target      <- as.numeric(get_arg("--q",         0.95))
stopifnot(map_model %in% c("standard", "hybrid", "synthetic"))

tilde_nu <- 3
## --s2-grid: explicit comma/space-separated s^2 values.  Overrides
## --s2-from/--s2-to/--s2-by when non-empty.  The default is a
## per-target dedicated cluster (3 points per N_target multiplier),
## union = 9 pts.  Centres are taken from the OBSERVED Pr-rule
## transitions in the previous run (not just a log-log fit):
##   * 1.00N (=100): 0.19, 0.20, 0.21  (linear interp of previous run
##                                      puts crossing at s^2 ~ 0.21;
##                                      0.16 was too far -- Ghat 0.46)
##   * 0.75N (=75):  0.22, 0.23, 0.24  (transition 0.94->0.97
##                                      between 0.22 and 0.23)
##   * 0.50N (=50):  0.28, 0.29, 0.30  (transition 0.89->0.98
##                                      between 0.28 and 0.30)
## (Pass --s2-grid "" to opt out and use the seq defaults.)
s2_grid_arg <- get_arg("--s2-grid",
                        "0.19,0.20,0.21,0.22,0.23,0.24,0.28,0.29,0.30")
if (nzchar(s2_grid_arg)) {
  S_grid <- sort(unique(as.numeric(strsplit(s2_grid_arg, "[, ]+")[[1]])))
} else {
  S_grid <- seq(s2_from, s2_to, by = s2_by)
}

## Reproducibility hygiene: clear THIS run's stale output(s) up front, so a
## crash (e.g. a multi-core PSOCK worker death -> "error reading from
## connection") can never leave a previous result on disk to be mistaken
## for fresh output.  Scoped to the file(s) this invocation writes (this
## map_model + RWD-size out_tag) -- sibling sizes (n200/n30) are untouched.
.resDir <- file.path(scriptDir, "res")
.stale  <- file.path(.resDir, c(
  sprintf("ess_calibration_%s%s%s%s.RData",            map_model, out_tag, cfg_tag, sc_suffix),
  sprintf("ess_calibration_%s%s%s%s_checkpoint.RData", map_model, out_tag, cfg_tag, sc_suffix)))
.stale  <- .stale[file.exists(.stale)]
if (length(.stale)) {
  cat("Clearing stale output before rerun:\n    ",
      paste(basename(.stale), collapse = "\n     "), "\n", sep = "")
  invisible(file.remove(.stale))
}

cat(sprintf("ess_calibration: B=%d  burn=%d  N_target=%g  q=%g  gMAP=%s  beta.prior=%s\n",
            B, burn, N_target, q_target, map_model, beta_prior_arg))
cat(sprintf("S_grid = [%s]  (K=%d points)\n",
            paste(S_grid, collapse=", "), length(S_grid)))

## --- canonical RWD ------------------------------------------------------

if (!file.exists(data_path))
  stop("--data path not found: '", data_path,
       "'  (omit --data to use the default: ", default_data, ")")
if (file.info(data_path)$isdir) {
  fs <- sort(list.files(data_path, pattern = "RData$", full.names = TRUE))
  stopifnot(length(fs) >= 1)
  canonical_file <- fs[1]
} else {
  canonical_file <- data_path
}
cat("Canonical RWD file:", basename(canonical_file), "\n")

d <- readRDS(canonical_file)
stopifnot("X" %in% names(d), "y" %in% names(d), "delta" %in% names(d),
          "D" %in% colnames(d$X), "Z" %in% colnames(d$X))

idx        <- d$X[, "D"] == 0 & d$X[, "Z"] == 0
Xcols      <- grep("^X[0-9]+$", colnames(d$X), value = TRUE)
X_RWD      <- as.matrix(d$X[idx, Xcols, drop = FALSE])
Yobs_RWD   <- as.numeric(d$y[idx])
delta_RWD  <- as.integer(d$delta[idx])
n_RWD      <- length(Yobs_RWD)
p          <- ncol(X_RWD)
log_Y      <- log(pmax(Yobs_RWD, .Machine$double.xmin))
sd_logY    <- sd(log_Y)
cat(sprintf("RWD subset (D==0 & Z==0): n=%d, p=%d, events=%d (%.0f%%), sd(log Y)=%.4f\n",
            n_RWD, p, sum(delta_RWD),
            100 * mean(delta_RWD), sd_logY))

## --- Stage 1: leaf summaries -------------------------------------------
## Run AFT-BART (heterbart + per-sweep rtnorm imputation) on
## (X_RWD, log_Y, delta_RWD); cabart_ess emits one row per
## (post-warmup draw b, tree h, leaf l): (b, h, leaf, ybar, n, sigma).

H        <- ntree_arg     # control-arm BART ntree (--ntree; default 50)
alpha    <- 0.95
beta     <- 2
k_tau    <- k_arg         # leaf-prior scale (--k; default 2): tau_theta = (max-min)/(2*k*sqrt(H))
numcut   <- 100L
nu       <- 3
sigquant <- 0.90

## AFT scale estimate via survreg on the RWD subset (matches mBART.R).
.aft_fit <- survreg(Surv(Yobs_RWD, delta_RWD) ~ ., data = data.frame(X_RWD),
                     dist = "lognormal")
sigest    <- .aft_fit$scale
lambda    <- (sigest^2 * qchisq(1 - sigquant, nu)) / nu
tau_theta <- (max(log_Y) - min(log_Y)) / (2 * k_tau * sqrt(H))
beta_prior <- if (identical(beta_prior_arg, "auto")) tau_theta else
                as.numeric(beta_prior_arg)
cat(sprintf("sigest(AFT scale)=%.3f  lambda=%.3f  tau_theta=%.3f  beta.prior=%.3f%s\n",
            sigest, lambda, tau_theta, beta_prior,
            if (identical(beta_prior_arg, "auto")) " (auto=tau_theta)" else ""))

source("/Users/oliviazhang/Desktop/bartModelMatrix.R")
temp    <- bartModelMatrix(X_RWD, numcut, usequants = FALSE, cont = FALSE,
                            xinfo = matrix(0, 0, 0), rm.const = TRUE)
Xt      <- t(temp$X)            # p x n
numcut  <- temp$numcut
xinfo   <- temp$xinfo
p_eff   <- nrow(Xt)
cat(sprintf("After bartModelMatrix: p_eff = %d (rm.const removed %d)\n",
            p_eff, p - p_eff))

t0 <- Sys.time()
cat("Compiling cabart_ess.cpp ... ")
sourceCpp(file.path(scriptDir, "cabart_ess.cpp"), verbose = FALSE)
cat(sprintf("done in %s\n", format(Sys.time() - t0)))

cat(sprintf("Running AFT-BART on RWD (H=%d, k=%g, B=%d, burn=%d) ...\n",
            H, k_tau, B, burn))
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

## sigma_ref: "auto" => use per-draw BART sigma_b (column sigma in
## leaf_table); numeric => override every row's sigma with that value.
use_auto_sigma <- identical(sigma_ref_arg, "auto")
sigma_ref      <- if (use_auto_sigma) NA_real_ else as.numeric(sigma_ref_arg)
cat(sprintf("sigma_ref = %s\n",
            if (use_auto_sigma) "auto: per-draw BART sigma_b (AFT log-scale)"
            else sprintf("%.4f (user-specified)", sigma_ref)))

## --- Stage 2: RBesT_lite gMAP + ELIR ESS -------------------------------

source(file.path(scriptDir, "RBesT_lite", "load.R"))

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

## One-row gMAP fit + automixfit + ess across the entire s^2 grid for a
## single (b, h, leaf) leaf summary.  Identical to the Gaussian script.
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
cat(sprintf("Per-worker progress lines (timestamp pid) follow.  Each is one b.\n\n"))
flush.console()

## PSOCK cluster, NOT mclapply: gMAP calls rstan::sampling, and
## rstan/RcppParallel/TBB thread state does not survive a fork() -- under
## mclapply the Stan sampler intermittently fails inside forked children,
## silently NA-ing a large fraction of draws.  PSOCK workers are fresh R
## processes; outfile="" forwards their progress to this terminal.
## CORES=1 runs sequentially in-process.
b_seq <- seq_len(B) - 1L
rbeslite_dir <- file.path(scriptDir, "RBesT_lite")
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
cat(sprintf("\nStage 2 done in %.1fs (%.1fm)  ~%.2fs / Stan fit\n",
            elapsed, elapsed / 60, elapsed * cores / total_cells))

## Pre-aggregate checkpoint: dump raw worker_results so a downstream
## aggregation crash doesn't cost the 4h Stage-2 budget.
resDir <- file.path(scriptDir, "res")
dir.create(resDir, showWarnings = FALSE)
.ckpt_path <- file.path(resDir,
        sprintf("ess_calibration_%s%s%s%s_checkpoint.RData", map_model, out_tag, cfg_tag, sc_suffix))
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
## Defensive: a worker may have died (Stan crash, OOM, etc.) and
## returned a try-error or NULL.  Drop those, warn, and aggregate over
## the survivors.

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
}
valid_results <- worker_results[ok_workers]
B_used <- length(valid_results)
stopifnot(B_used > 0)

ESSbar <- do.call(rbind, lapply(valid_results, `[[`, "ESSbar_b"))  # B_used x K (rbind is robust to K=1; t(sapply(...)) transposes to 1xB there)

## Derive c_s2 / Ghat directly from ESSbar with na.rm so a draw whose
## ESSbar_b is NA in one s^2 column does not poison the whole column.
below_mat <- ESSbar <= N_target
c_s2      <- colSums(below_mat, na.rm = TRUE)
n_valid   <- colSums(!is.na(ESSbar))
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
## Plotting is delegated to ess_plot.R (or a future ess_plot_surv.R).

out <- list(
  canonical_file = canonical_file,
  scenario       = sc_tag,
  n_RWD          = n_RWD,
  p              = p_eff,
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
  tilde_nu       = tilde_nu,
  map_model   = map_model,
  beta_prior     = beta_prior,
  sigma_ref      = sigma_ref,
  outcome        = "surv"
)
saveRDS(out, file = file.path(resDir,
        sprintf("ess_calibration_%s%s%s%s.RData", map_model, out_tag, cfg_tag, sc_suffix)))

cat("Wrote", file.path(resDir,
       sprintf("ess_calibration_%s%s%s%s.RData", map_model, out_tag, cfg_tag, sc_suffix)), "\n")
cat("Run ess_plot.R --outcome surv to generate the figure.\n")
