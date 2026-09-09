## ess_calibration_map.R -- Prior ESS calibration for the *vanilla*
## MAP prior (as in MAP.R).  No BART, no trees: the canonical RWD's
## D == 0 & Z == 0 subset is summarised as a single (ybar, n, SE) row
## and the MAP prior is fit with RBesT_lite::gMAP exactly as MAP.R
## does it.
##
## Companion to ess_calibration.R (MAP-BART): the two scripts share CLI
## flags, section order, and the same Stage-2 gMAP + automixfit + ELIR
## ESS core; this script differs only in Stage 1, where the "leaf
## table" collapses to a single row (the whole RWD subset).  Run both
## on the same RWD to compare vanilla-MAP vs MAP-BART calibration on
## identical settings.
##
## Because the leaf table is deterministic given the RWD, ESS(s^2) is
## one number per s^2 (modulo small MCMC noise inside gMAP), so the
## Pr-rule used by ess_calibration.R collapses to a single threshold:
##   hat_s2_star = min{ s^2 in grid : ESS(s^2) <= N_target }.
##
## Usage:
##   Rscript ess_calibration_map.R [--data <RWD file or dir>]
##                                  [--s2-from 0.05] [--s2-to 0.45] [--s2-by 0.05]
##                                  [--s2-grid "0.05,0.1,..."]   # uneven; overrides from/to/by
##                                  [--Ntarget 100] [--seed 6]
##                                  [--map-model standard]
##                                  [--beta-prior 10] [--sigma-ref auto]
##
## CLI is a strict subset of ess_calibration.R: the flags above are
## shared (same order, possibly different defaults).  The BART-only
## flags (--B, --burn, --q) do not appear here -- there are no MCMC
## draws to average over and no Pr-rule to apply.
##
## --map-model selects which RBesT_lite/stan/gMAP_*.stan is fit:
##   standard (default)  -- the default MAP: theta_h = mu + tau * xi_h,
##                          mu ~ N(0, beta.prior^2) (matches MAP.R)
##   hybrid              -- two-arm hybrid design: theta_h ~ N(0, tau^2),
##                          mu fixed at 0 (zero-centred MAP model, §2.3;
##                          beta.prior is ignored)
##   synthetic           -- single-arm synthetic control (single anchor)
##
## --sigma-ref controls the sigma passed to ess(MAP_mixture, sigma=.):
##   auto (default)     -- sd(y) of the RWD subset, so ESS reads as
##                          "the prior carries the information of [ESS]
##                          outcome-scale observations"
##   <numeric>          -- fixed user-specified value
##
## Defaults: canonical RWD = first *.RData file (alphabetical) under
## /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/data_v2/, D == 0 & Z == 0
## subset.

scriptDir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)),
  error = function(e) "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/ess_local"
)
if (!nzchar(scriptDir) || is.na(scriptDir))
  scriptDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/ess_local"

## --- CLI ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  k <- grep(paste0("^", flag, "$"), args)
  if (length(k) == 0) return(default)
  args[k + 1]
}

default_data  <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/data_v2"
data_path     <- get_arg("--data",        default_data)
s2_from       <- as.numeric(get_arg("--s2-from",   0.05))
s2_to         <- as.numeric(get_arg("--s2-to",     0.45))
s2_by         <- as.numeric(get_arg("--s2-by",     0.05))
N_target      <- as.numeric(get_arg("--Ntarget",   100))
seed          <- as.integer(get_arg("--seed",      6L))
map_model  <- get_arg("--map-model", "standard")
beta_prior    <- as.numeric(get_arg("--beta-prior", 10))
sigma_ref_arg <- get_arg("--sigma-ref", "auto")   # "auto" = sd(y) of RWD subset
stopifnot(map_model %in% c("standard", "hybrid", "synthetic"))

tilde_nu <- 3
## --s2-grid: explicit comma/space-separated s^2 values; overrides
## --s2-from/--s2-to/--s2-by when given.
s2_grid_arg <- get_arg("--s2-grid", "")
if (nzchar(s2_grid_arg)) {
  S_grid <- sort(unique(as.numeric(strsplit(s2_grid_arg, "[, ]+")[[1]])))
} else {
  S_grid <- seq(s2_from, s2_to, by = s2_by)
}

## Reproducibility hygiene: clear THIS run's stale output up front, so a
## crash can never leave a previous result on disk to be mistaken for fresh
## output.  Scoped to the file this invocation writes.
.resDir <- file.path(scriptDir, "res")
.stale  <- file.path(.resDir, sprintf("ess_calibration_map_%s.RData", map_model))
.stale  <- .stale[file.exists(.stale)]
if (length(.stale)) {
  cat("Clearing stale output before rerun:\n     ", basename(.stale), "\n", sep = "")
  invisible(file.remove(.stale))
}

cat(sprintf("ess_calibration_map: N_target=%g  gMAP=%s  beta.prior=%g\n",
            N_target, map_model, beta_prior))
cat(sprintf("S_grid = [%s]  (K=%d points)\n",
            paste(S_grid, collapse=", "), length(S_grid)))

## --- canonical RWD -----------------------------------------------------

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

d   <- readRDS(canonical_file)
idx <- d$X[, "D"] == 0 & d$X[, "Z"] == 0
y   <- d$y[idx]
n_RWD <- length(y)
ybar  <- mean(y)
sd_y  <- sd(y)
SE    <- sd_y / sqrt(n_RWD)
cat(sprintf("RWD subset (D==0 & Z==0): n=%d, ybar=%.4f, sd(y)=%.4f, SE=%.4f\n",
            n_RWD, ybar, sd_y, SE))

## --- Stage 1: leaf summaries -------------------------------------------
## Vanilla MAP: one "leaf" = the whole RWD subset, so the leaf table has
## exactly one row.  Column names match the BART version's leaf_table.

leaf_table <- data.frame(b = 0L, h = 0L, leaf = 0L,
                          ybar = ybar, n = n_RWD, sigma = sd_y)

## Default sigma_ref: sd(y) of the RWD subset, i.e., the natural
## noise scale of one observation.  ESS then reads as "the prior
## carries the information of [ESS] outcome-scale observations,"
## directly comparable to n_RWD.
sigma_ref <- if (identical(sigma_ref_arg, "auto")) sd_y else
              as.numeric(sigma_ref_arg)
cat(sprintf("sigma_ref = %.4f  (%s)\n",
            sigma_ref,
            if (identical(sigma_ref_arg, "auto")) "auto: sd(y) of RWD subset"
            else "user-specified"))

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

one_row <- data.frame(y.bar = ybar, n = n_RWD, y.se = SE,
                       study = factor(1, levels = "1"))

set.seed(seed)
K   <- length(S_grid)
ESS <- numeric(K)

t0 <- Sys.time()
for (k in seq_along(S_grid)) {
  s2 <- S_grid[k]
  ts <- Sys.time()
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
  ESS[k]      <- ess(MAP_mixture, sigma = sigma_ref)
  message(sprintf("[%s] s^2=%-7.4g  ESS=%-8.2f  (%.2fs)",
                  format(Sys.time(), "%H:%M:%S"),
                  s2, ESS[k], as.numeric(Sys.time() - ts)))
}
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("Total %.1fs across K=%d s^2 values\n", elapsed, K))

## --- Aggregate + selection rule ----------------------------------------

satisfies   <- which(ESS <= N_target)
hat_s2_star <- if (length(satisfies)) S_grid[min(satisfies)] else NA_real_

cat("\n--- Vanilla-MAP ESS calibration results ---\n")
print(data.frame(s2 = S_grid, ESS = round(ESS, 2),
                  below_target = ESS <= N_target))
if (is.na(hat_s2_star)) {
  cat(sprintf("\nNo s^2 in grid satisfies ESS <= %g; extend grid upward.\n",
              N_target))
} else {
  cat(sprintf("\nCalibrated s^2_star = %.4g  (smallest with ESS <= %g)\n",
              hat_s2_star, N_target))
}

## --- Save ---------------------------------------------------------------
## Plotting is delegated to ess_plot.R, which loads this .RData and the
## BART .RData and produces a single merged figure.

resDir <- file.path(scriptDir, "res")
dir.create(resDir, showWarnings = FALSE)

out <- list(
  canonical_file = canonical_file,
  n_RWD          = n_RWD,
  ybar           = ybar,
  SE             = SE,
  S_grid         = S_grid,
  ESS            = ESS,
  leaf_table     = leaf_table,
  N_target       = N_target,
  hat_s2_star    = hat_s2_star,
  tilde_nu       = tilde_nu,
  map_model   = map_model,
  beta_prior     = beta_prior,
  sigma_ref      = sigma_ref
)
saveRDS(out, file = file.path(resDir,
        sprintf("ess_calibration_map_%s.RData", map_model)))

cat("Wrote", file.path(resDir,
       sprintf("ess_calibration_map_%s.RData", map_model)), "\n")
cat("Run ess_plot.R to generate the merged figure.\n")
