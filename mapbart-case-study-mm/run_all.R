## mapbart-case-study-mm/run_all.R
##
## Runs the four analysis scripts -- km_analysis.R, AFTv2_analysis.R,
## HierAFT_analysis.R, mBART_analysis.R -- for a single OUTCOME (PFS
## default, or OS via OUTCOME=OS), then loads their saved result files
## and prints a single results table (delta_hat + 95% CI per method /
## per N_target).
##
## Each script reads OUTCOME via Sys.getenv("OUTCOME", unset = "PFS"),
## so we invoke each in a fresh Rscript subprocess.  Using subprocesses
## avoids state leakage (rstan/RcppParallel thread state, compiled C++,
## large objects) between scripts.
##
## Output files read for the table (OUTCOME=PFS, n283 shown):
##   KM_results_PFS_n283.RData
##   AFTv2_results_PFS_n283.RData
##   HierAFT_results_PFS_n283.RData
##   mBART_results_PFS_n283_N<target>.RData    (one per N_target multiplier)
##
## OUTCOME=OS uses the matching OS calibration
## (ess_calibration_ucmm_os_*_n<N>.RData); mBART's discrepancy-prior rate
## is set at runtime from the calibrated s^2 (cmbart setrate), so no C++
## recompile is needed to switch outcomes.
##
## Usage:
##   Rscript run_all.R                              # default (PFS)
##   OUTCOME=OS Rscript run_all.R                   # overall survival
##   SCRIPTS="mBART_analysis.R" Rscript run_all.R   # rerun a subset
##   MBART_CONFIGS="default" Rscript run_all.R      # pick control-arm BART
##                                                  # config(s): names
##                                                  # shrunk / mid / default,
##                                                  # or ntree:k pairs;
##                                                  # comma-separated, e.g.
##                                                  # "shrunk,default"
##   SKIP_RUN=1 Rscript run_all.R                   # just rebuild the
##                                                  # table from existing
##                                                  # res/ files

rm(list = ls())

mainDir <- "/Users/oliviazhang/Desktop"
projDir <- file.path(mainDir, "mapbart-case-study-mm")
res_dir <- file.path(projDir, "res")

## OUTCOME: PFS (default) or OS.  Propagated to every subprocess and used
## to tag the result files it reads back.
OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))

## Merged cohort selected by FILE NAME. Override per run, e.g. the secondary
## (KRd/Rd) cohort:
##   MERGED_FILE=/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n78.RData Rscript run_all.R
## Default = primary cohort: 30 EloKRD + 253 UCMM (all regimens, E-Rd excluded).
MERGED_FILE <- Sys.getenv("MERGED_FILE",
                          unset = file.path(projDir, "data_cleaned/merged_elokrd_ucmm_n283.RData"))
## Allow a bare filename (as the README and analysis scripts do): resolve it
## against data_cleaned/ so run_all.R's own load() below finds it too.
if (!file.exists(MERGED_FILE) &&
    file.exists(file.path(projDir, "data_cleaned", basename(MERGED_FILE))))
  MERGED_FILE <- file.path(projDir, "data_cleaned", basename(MERGED_FILE))
if (!file.exists(MERGED_FILE))
  stop("Merged file not found: ", MERGED_FILE, "  (run data_merge.R first)")
data_tag <- regmatches(basename(MERGED_FILE), regexpr("n[0-9]+", basename(MERGED_FILE)))
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(MERGED_FILE))

## --- CLI / env-var configuration ----------------------------------------

default_scripts <- c("km_analysis.R",
                      "AFTv2_analysis.R",
                      "HierAFT_analysis.R",
                      "BARTv2_analysis.R",
                      "mBART_analysis.R")
scripts <- strsplit(Sys.getenv("SCRIPTS",
                                 unset = paste(default_scripts, collapse=",")),
                     "[, ]+")[[1]]
scripts <- scripts[nzchar(scripts)]
## SKIP_RUN accepts the usual truthy strings (1, TRUE, true, T, yes);
## anything else (or unset) is false.  Avoids as.logical("1") -> NA.
skip_run <- Sys.getenv("SKIP_RUN", unset = "") %in%
              c("1", "TRUE", "true", "T", "yes", "Y", "y")

## --- Control-arm BART config selection for mBART_analysis.R -------------
## MBART_CONFIGS chooses which control-arm (ntree, k) config(s) mBART runs,
## as comma-separated explicit ntree:k pairs:
##   MBART_CONFIGS="50:2"           # one config
##   MBART_CONFIGS="5:2,10:2,5:4"   # several configs
## Unset -> the default sweep below.  Passed to mBART via the MBART_CONFIGS
## env var; the same pairs drive which config-tagged result files get
## aggregated below.
mbart_default  <- "5:2,10:2,5:4"   # default sweep when MBART_CONFIGS is unset
.cfg_tokens <- strsplit(Sys.getenv("MBART_CONFIGS", unset = mbart_default),
                        "[, ]+")[[1]]
.cfg_tokens <- .cfg_tokens[nzchar(.cfg_tokens)]
.resolve_cfg <- function(tok) {
  if (grepl("^[0-9]+:[0-9.]+$", tok)) return(tok)
  stop("Bad MBART_CONFIGS token '", tok, "'  (expected ntree:k pairs, e.g. 50:2)")
}
mbart_pairs <- unique(vapply(.cfg_tokens, .resolve_cfg, character(1)))
mbart_configs_wire <- paste(mbart_pairs, collapse = ",")

## --- MAP-BART s^2 at each N_target multiplier: from the ess calibration ------
## Read the same calibration mBART_analysis.R uses (OUTCOME + data_tag) and, for
## each candidate multiplier m, apply the ess Pr-rule to extract the calibrated
## leaf-discrepancy s^2:
##     s^2(m) = smallest s^2 on S_grid with  Pr{ ESSbar <= m x N_target } >= q.
## Keep only multipliers whose Ghat(m x N_target) crosses q at an INTERIOR grid
## point -- it reaches q AND is not already >= q at the smallest s^2 (i.e. not
## saturated, so .pick returns a finite s^2).  This is the SAME rule used by
## mapbart-case-study-mm/mBART_analysis.R's .pick and by the mapbart-sim-*/run_all.R +
## mapbart-sim-*/mBART.R pair (run_all picks the multipliers; mBART picks the
## s^2) -- here we do both so run_all reports the s^2 it is calibrating to.
## mbart_mult is passed to mBART_analysis.R via MBART_NMULT so the analysis and
## this aggregation agree on which targets exist; mBART recomputes s^2 from the
## same rule.
mbart_mult_pool <- c(3, 2, 1, 0.75, 0.50)
.cal_dir <- file.path(projDir, "ess_local", "res")
# Config-matched calibration file: REQUIRE the exact _nt<ntree>_k<k> tag, same as
# mBART_analysis.R / mapbart-sim-survival-single-arm -- one calibration per control config.
.cal_file_for <- function(nt, k) {
  pat <- sprintf("^ess_calibration_ucmm_%s_.*_%s_nt%d_k%g\\.RData$", tolower(OUTCOME), data_tag, nt, k)
  fs  <- list.files(.cal_dir, pattern = pat, full.names = TRUE)
  fs  <- fs[!grepl("_checkpoint\\.RData$", fs)]
  if (!length(fs)) return(NA_character_)
  fs[order(file.info(fs)$mtime, decreasing = TRUE)][1]
}
# Multipliers a calibration supports (interior q-crossing only).  Same rule as
# mBART_analysis.R's .pick and mapbart-sim-survival-single-arm/run_all.R's mbart_mult_from.
.mult_from <- function(cal) {
  if (is.null(cal)) return(numeric(0))
  av <- Filter(function(m) {
    Gh <- colMeans(cal$ESSbar <= m * cal$N_target, na.rm = TRUE)
    any(Gh >= cal$q) && Gh[1] < cal$q
  }, mbart_mult_pool)
  sort(unlist(av), decreasing = TRUE)
}
# Each control config reads its OWN calibration; MAP-BART runs a config only at
# the N that config supports (mBART_analysis.R does the per-config s^2 pick).
# The UNION of supported target Ns across the AVAILABLE configs is what AFTv2
# sweeps (so AFTv2 runs at every N any config supports) and what labels the
# plot.  Mirrors mapbart-sim-survival-single-arm/run_all.R (ess_mults / ess_target_Ns).
.cfg_list <- lapply(strsplit(mbart_pairs, ":"),
                    function(v) list(ntree = as.integer(v[1]), k = as.numeric(v[2])))
.all_mult <- numeric(0); .all_tn <- integer(0)
for (cf in .cfg_list) {
  f <- .cal_file_for(cf$ntree, cf$k)
  if (is.na(f)) {
    cat(sprintf("mapbart-case-study-mm/run_all.R: no %s/%s calibration for config nt%d_k%g -- MAP-BART will skip it\n",
                tolower(OUTCOME), data_tag, cf$ntree, cf$k)); next
  }
  cal <- readRDS(f); ms <- .mult_from(cal)
  tns <- as.integer(round(ms * cal$N_target))
  .all_mult <- c(.all_mult, ms); .all_tn <- c(.all_tn, tns)
  cat(sprintf("mapbart-case-study-mm/run_all.R: config nt%d_k%g <- %s: supported N = %s\n",
              cf$ntree, cf$k, basename(f),
              if (length(tns)) paste(tns, collapse = ",") else "(none)"))
}
mbart_mult      <- sort(unique(.all_mult), decreasing = TRUE)
mbart_target_Ns <- sort(unique(.all_tn),   decreasing = TRUE)
if (!length(mbart_mult)) {
  mbart_mult      <- c(1.00, 0.75, 0.50)
  mbart_target_Ns <- integer(0)
  warning("mapbart-case-study-mm/run_all.R: no config-matched ess calibration found; mBART multipliers fall back to 1,0.75,0.5")
}
mbart_nmult_wire <- paste(mbart_mult, collapse = ",")

## --- AFTv2 power-prior weights (rwd_w): matched effective external-control N --
## Mirrors mapbart-sim-survival-single-arm/run_all.R's aftv2_rwd_w_vals: rwd_w = 1 (full
## borrowing) plus N_target / n_rwd for each finite MAP-BART target N, so AFTv2
## is compared at a matched effective external-control sample size.  n_rwd is the
## ACTUAL UCMM (external-control) count for this cohort (read from MERGED_FILE),
## not a hard-coded nominal, since mapbart-case-study-mm's design is external-control-only.
.merged_env <- new.env()
load(MERGED_FILE, envir = .merged_env)
n_rwd_nominal <- sum(.merged_env$merged$trt == 0)
aftv2_rwd_w_vals <- sort(unique(c(1, if (length(mbart_target_Ns))
                                       round(mbart_target_Ns / n_rwd_nominal, 4))),
                          decreasing = TRUE)
aftv2_rwd_w_wire <- paste(aftv2_rwd_w_vals, collapse = ",")
cat(sprintf("mapbart-case-study-mm/run_all.R: AFTv2 rwd_w values = %s  (n_UCMM = %d)\n",
            aftv2_rwd_w_wire, n_rwd_nominal))

## --- Prior sweep for HierAFT_analysis.R ---------------------------------
## HIERAFT_CONFIGS chooses which discrepancy-variance prior(s) HierAFT runs,
## as comma-separated values for IG(3/2, 3*prior/2):
##   HIERAFT_CONFIGS="0.25"          # one prior
##   HIERAFT_CONFIGS="0.1,0.25,0.5"  # a sweep
## Unset -> the default sweep below.  Passed to HierAFT via the HIERAFT_CONFIGS
## env var; the same values drive which prior-tagged result files get read.
hieraft_default <- "0.05,0.5"   # default sweep when HIERAFT_CONFIGS is unset
                                # (mirrors mapbart-sim-survival-single-arm/run_all.R lmv4_prior_vals = c(0.05, 0.5))
.haft_tokens   <- strsplit(Sys.getenv("HIERAFT_CONFIGS", unset = hieraft_default), "[, ]+")[[1]]
hieraft_priors <- as.numeric(.haft_tokens[nzchar(.haft_tokens)])
hieraft_priors <- hieraft_priors[!is.na(hieraft_priors)]
hieraft_configs_wire <- paste(hieraft_priors, collapse = ",")

cat(sprintf("mapbart-case-study-mm/run_all.R: OUTCOME = %s   MERGED_FILE = %s   scripts = %s\n   mBART configs = %s   HierAFT priors = %s%s\n",
            OUTCOME, basename(MERGED_FILE), paste(scripts, collapse=", "),
            sprintf("%s   mBART N-mult = %s", mbart_configs_wire, mbart_nmult_wire), hieraft_configs_wire,
            if (skip_run) "   [SKIP_RUN: only rebuild table]" else ""))

## --- run each script ----------------------------------------------------

run_one <- function(script) {
  banner <- strrep("=", 64)
  cat(sprintf("\n%s\n=== %s   (OUTCOME=%s) ===\n%s\n",
              banner, script, OUTCOME, banner))
  flush.console()
  t0 <- Sys.time()
  status <- system2("Rscript",
                     args = shQuote(file.path(projDir, script)),
                     env  = c(sprintf("OUTCOME=%s", OUTCOME),
                              sprintf("MERGED_FILE=%s", MERGED_FILE),
                              sprintf("MBART_CONFIGS=%s", mbart_configs_wire),
                              sprintf("MBART_NMULT=%s", mbart_nmult_wire),
                              sprintf("AFTV2_RWD_W=%s", aftv2_rwd_w_wire),
                              sprintf("HIERAFT_CONFIGS=%s", hieraft_configs_wire)))
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf(">>> %s [%s]: exit=%d  elapsed=%.1fs\n",
              script, OUTCOME, status, dt))
  status
}

if (!skip_run) {
  for (s in scripts) run_one(s)
}

## --- aggregate the results into one table ------------------------------
## Each script's result list has $delta_hat (point estimate) and
## $ci_95 (length-2 vector: lower, upper).  For mBART there is one file
## per N_target multiplier; we read all matching files.

extract_row <- function(method, path) {
  if (!file.exists(path)) {
    return(data.frame(method = method, target_N = NA_integer_, prior = NA_real_,
                       delta_hat = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                       delta_diff = NA_real_, diff_lo = NA_real_, diff_hi = NA_real_,
                       rmst_trt = NA_real_,  rmst_trt_lo = NA_real_,  rmst_trt_hi = NA_real_,
                       rmst_ctrl = NA_real_, rmst_ctrl_lo = NA_real_, rmst_ctrl_hi = NA_real_,
                       sigma_trt = NA_real_, sigma_trt_lo = NA_real_, sigma_trt_hi = NA_real_,
                       sigma_ctrl = NA_real_, sigma_ctrl_lo = NA_real_, sigma_ctrl_hi = NA_real_,
                       file = basename(path), stringsAsFactors = FALSE))
  }
  r <- readRDS(path)
  dh <- if (!is.null(r$delta_hat)) r$delta_hat else NA_real_
  ci <- if (!is.null(r$ci_95) && length(r$ci_95) == 2) r$ci_95
        else c(NA_real_, NA_real_)
  tN <- if (!is.null(r$settings) && !is.null(r$settings$target_N))
          as.integer(r$settings$target_N)
        else NA_integer_
  ## prior s^2 (mBART discrepancy prior); used to label the fixed-s^2 run.
  pr  <- if (!is.null(r$settings) && !is.null(r$settings$prior))
           as.numeric(r$settings$prior)
         else NA_real_
  dd  <- if (!is.null(r$delta_diff)) r$delta_diff else NA_real_
  cid <- if (!is.null(r$ci_diff_95) && length(r$ci_diff_95) == 2) r$ci_diff_95
         else c(NA_real_, NA_real_)
  ## arm-specific estimates (population RMST + residual sigma per arm; est + 95% CI)
  gv <- function(v, k) if (!is.null(v)) unname(v[k]) else NA_real_
  data.frame(method = method, target_N = tN, prior = pr,
             delta_hat = unname(dh),
             ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
             delta_diff = unname(dd), diff_lo = unname(cid[1]), diff_hi = unname(cid[2]),
             rmst_trt  = gv(r$rmst_trt_est,"est"),  rmst_trt_lo = gv(r$rmst_trt_est,"lo"),  rmst_trt_hi = gv(r$rmst_trt_est,"hi"),
             rmst_ctrl = gv(r$rmst_ctrl_est,"est"), rmst_ctrl_lo = gv(r$rmst_ctrl_est,"lo"), rmst_ctrl_hi = gv(r$rmst_ctrl_est,"hi"),
             sigma_trt  = gv(r$sigma_trt_est,"est"),  sigma_trt_lo = gv(r$sigma_trt_est,"lo"),  sigma_trt_hi = gv(r$sigma_trt_est,"hi"),
             sigma_ctrl = gv(r$sigma_ctrl_est,"est"), sigma_ctrl_lo = gv(r$sigma_ctrl_est,"lo"), sigma_ctrl_hi = gv(r$sigma_ctrl_est,"hi"),
             file = basename(path), stringsAsFactors = FALSE)
}

table_rows <- list()
## Row 1: raw KM RMST(5yr) ratio (empirical anchor, no borrowing) -- produced
## by km_analysis.R, read here like any other method's result file.
table_rows[[1]] <- extract_row("KM",
                                file.path(res_dir, sprintf("KM_results_%s_%s.RData", OUTCOME, data_tag)))
## AFTv2: one result file per power-prior weight rwd_w, tagged _w<rwd_w>.
## extract_row pulls each run's matched target_N from settings for the N column.
for (rw in aftv2_rwd_w_vals) {
  table_rows[[length(table_rows) + 1]] <-
    extract_row(sprintf("AFTv2(w=%g)", rw),
                file.path(res_dir, sprintf("AFTv2_results_%s_%s_w%g.RData", OUTCOME, data_tag, rw)))
}
table_rows[[length(table_rows) + 1]] <-
  extract_row("BARTv2",
              file.path(res_dir, sprintf("BARTv2_results_%s_%s.RData", OUTCOME, data_tag)))
## HierAFT: one result file per prior config, tagged _prior<prior>.
for (pr in hieraft_priors) {
  table_rows[[length(table_rows) + 1]] <-
    extract_row(sprintf("HierAFT(prior=%g)", pr),
                file.path(res_dir, sprintf("HierAFT_results_%s_%s_prior%g.RData", OUTCOME, data_tag, pr)))
}

## mBART: one set of result files per selected control-arm config, tagged
## _nt<ntree>_k<k>.  Label each row with its config so configs are
## distinguishable in the table.
for (pair in mbart_pairs) {
  v  <- as.numeric(strsplit(pair, ":")[[1]])
  nt <- as.integer(v[1]); kk <- v[2]
  pat <- sprintf("^mBART_results_%s_%s_nt%d_k%g_N[0-9]+\\.RData$",
                 OUTCOME, data_tag, nt, kk)
  files <- list.files(res_dir, pattern = pat, full.names = TRUE)
  files <- files[order(file.info(files)$mtime)]
  for (f in files)
    table_rows[[length(table_rows) + 1]] <-
      extract_row(sprintf("mBART(nt%d,k%g)", nt, kk), f)
}

results_table <- do.call(rbind, table_rows)
## Pull target_N out of the mBART filename when settings$target_N is
## absent (older saves).
.from_name <- function(f) {
  m <- regmatches(f, regexpr("_N[0-9]+", f))
  if (length(m)) as.integer(sub("_N", "", m)) else NA_integer_
}
need <- is.na(results_table$target_N) & grepl("^mBART", results_table$method)
if (any(need))
  results_table$target_N[need] <- vapply(results_table$file[need], .from_name, integer(1))

## Widen print width so the formatted "est [lo, hi]" tables print one row
## per method instead of wrapping the trailing columns into a second block.
options(width = 200)

## N label: no-borrowing references show "-"; the fixed-s^2 mBART run
## (target_N sentinel 0) shows its s^2 instead of a misleading "0".
.Nlab <- function(tN, pr) ifelse(is.na(tN), "-",
                                 ifelse(tN == 0L, sprintf("s2=%g", pr), as.character(tN)))
fmt_ci <- function(est, lo, hi)
  ifelse(is.na(est), "NA", sprintf("%.3f [%.3f, %.3f]", est, lo, hi))

## --- (1) Treatment-effect estimates: RMST(5yr) ratio + difference -------
cat("\n", strrep("=", 64),
    "\n=== Treatment effect (OUTCOME=", OUTCOME, ", data=", data_tag, ") ===\n",
    strrep("=", 64), "\n", sep = "")
disp_eff <- with(results_table, data.frame(
  method     = method,
  N          = .Nlab(target_N, prior),
  rmst_ratio = fmt_ci(delta_hat,  ci_lo,   ci_hi),
  rmst_diff  = fmt_ci(delta_diff, diff_lo, diff_hi),
  stringsAsFactors = FALSE, check.names = FALSE))
print(disp_eff, row.names = FALSE)

## --- (2) Arm-specific estimates: population RMST + residual sigma per arm
cat("\n", strrep("-", 64),
    "\n=== Arm-specific estimates (population RMST and residual sigma) ===\n",
    strrep("-", 64), "\n", sep = "")
disp_arm <- with(results_table, data.frame(
  method     = method,
  N          = .Nlab(target_N, prior),
  rmst_trt   = fmt_ci(rmst_trt,   rmst_trt_lo,   rmst_trt_hi),
  rmst_ctrl  = fmt_ci(rmst_ctrl,  rmst_ctrl_lo,  rmst_ctrl_hi),
  sigma_trt  = fmt_ci(sigma_trt,  sigma_trt_lo,  sigma_trt_hi),
  sigma_ctrl = fmt_ci(sigma_ctrl, sigma_ctrl_lo, sigma_ctrl_hi),
  stringsAsFactors = FALSE, check.names = FALSE))
print(disp_arm, row.names = FALSE)
cat("\nrmst_ratio/diff = RMST(5yr) trt-vs-ctrl; rmst_*/sigma_* per arm; cells = est [95% CI].",
    "\nN = mBART target_N (s2=... = fixed-prior run); '-' = no-borrowing reference.\n")

saveRDS(results_table,
        file = file.path(res_dir, sprintf("run_all_summary_%s_%s.RData", OUTCOME, data_tag)))
cat(sprintf("\nSaved aggregated table to: %s\n",
            file.path(res_dir, sprintf("run_all_summary_%s_%s.RData", OUTCOME, data_tag))))
