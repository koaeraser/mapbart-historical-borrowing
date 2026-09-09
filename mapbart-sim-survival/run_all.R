# Master script: full pipeline.
#
#   1. Generate data    via  data_gen_p<p_obs>_v2.R
#   2. Plot balance     via  plot_balance.R
#   3. Run analyses     via  AFTv1, AFTv2, AFTv3, HierAFT, BARTv1, BARTv2, BARTv3, mBART
#   4. Plot results     via  plot.R
#
# The master reads each file's text in memory, regex-replaces the matching
# header assignments with the current overrides, and evaluates the modified
# text in an isolated environment. The files on disk are never modified.
#
# `p_obs` is auto-inferred inside each analysis file from the data's $X
# columns, so we deliberately do NOT pass it as an override there. It is
# passed to data_gen and the plot scripts (which use it directly).

rm(list = ls())

mainDir <- "/Users/oliviazhang/Desktop/"
script_dir <- file.path(mainDir, "mapbart-sim-survival")

# ---- Pipeline flags ----
run_data_gen     <- TRUE
run_plot_balance <- TRUE
run_analysis     <- TRUE
run_plot         <- TRUE

# Run analysis methods in parallel within each scenario (mclapply, fork-based).
# Set FALSE to run sequentially; falls back to sequential on Windows.
parallel_methods <- TRUE
use_parallel     <- isTRUE(parallel_methods) && .Platform$OS.type == "unix"  # parallelize methods (unix only)
n_cores          <- NULL   # NULL => use length(files); otherwise set an integer

# ---- Scenario sweep grid ----
p_obs        <- 10
sc_values    <- c(1, 2, 3)
subsc_values <- c("E")
hypo_values  <- c("alternative")  # add "null" first to calibrate thresholds
cor_for_sc   <- function(sc) if (sc == 3) c(-0.5, 0, 0.5) else 1

# ---- Analysis files ----
files <- c("AFTv1.R", "AFTv2.R", "AFTv3.R", "BARTv1.R", "BARTv2.R", "BARTv3.R", "HierAFT.R", "mBART.R")
# Files in `serial_files` run sequentially AFTER the parallel batch (use for
# methods whose native libraries aren't fork-safe). rstan breaks when several
# Stan models are mcparallel-forked at once -> "invalid connection". So the
# Stan-based AFT methods run serially; only the fork-safe ones (BART C++) go
# in the parallel batch.
serial_files <- c("AFTv1.R", "AFTv2.R", "AFTv3.R", "HierAFT.R")

# ---- Seeds ----
# seed_rct      : data_gen RCT iteration seed-array generator (rewrites set.seed(<literal>))
# seed_rwd      : data_gen RWD-pool seed (overrides `seed_rwd <- ...` in data_gen)
# seed_method : analysis-script seed-array generator (rewrites set.seed(<literal>))
seed_rct      <- 123
seed_rwd      <- 456
seed_method <- 789

# ---- Shared overrides ----
# Applied to all sourced files (no-op on files that don't have these vars).
common_overrides <- list(
  seed_method = seed_method,  # seed-array generator for analysis scripts
  rmst_control_sigma = "own"          # control-arm RMST sigma: "trt"|"own"|"adaptive" (robust across n)
)
# Per-file overrides, layered on top of common_overrides and the scenario
# overrides inside `run_one()`.
#
# mBART.R no longer needs `prior_vals` here: it loads
# ess_local/res/ess_calibration_surv_<version>.RData itself and iterates
# over the N_target multipliers (1.00, 0.75, 0.50 by default),
# labelling output files by the integer target N rather than the raw
# s^2 value.  HierAFT.R is the parametric AFT analog (no leaf prior to
# calibrate), so it keeps a fixed prior_vals override.
# ---- MAP-BART N targets: extracted from the ess calibration results ----
# From a candidate multiplier pool, keep each m whose Ghat(m x N_target) crosses
# q at an INTERIOR grid point (reaches q AND not already >= q at the smallest
# s^2 -> not saturated).  This auto-detects which targets the calibration
# supports instead of hand-coding multipliers.  mbart_mult is injected into
# mBART.R (N_target_multipliers -> which targets it runs) and the plot (labels).
mbart_mult_pool <- c(3, 2, 1, 0.75, 0.50)
.cal_files <- list.files(file.path(script_dir, "ess_local", "res"),
                         pattern = "^ess_calibration_[^/]*\\.RData$", full.names = TRUE)
.cal_files <- .cal_files[!grepl("_checkpoint\\.RData$", .cal_files)]
if (length(.cal_files)) {
  .cal <- readRDS(.cal_files[order(file.info(.cal_files)$mtime, decreasing = TRUE)][1])
  .avail <- Filter(function(m) {
    Gh <- colMeans(.cal$ESSbar <= m * .cal$N_target, na.rm = TRUE)
    any(Gh >= .cal$q) && Gh[1] < .cal$q
  }, mbart_mult_pool)
  mbart_mult      <- sort(unlist(.avail), decreasing = TRUE)
  mbart_target_Ns <- c(Inf, as.integer(round(mbart_mult * .cal$N_target)))
  cat(sprintf("run_all.R: MAP-BART multipliers from ess results = {%s} -> targets %s\n",
              paste(mbart_mult, collapse = ", "), paste(mbart_target_Ns, collapse = ", ")))
} else {
  mbart_mult      <- c(1.00, 0.75, 0.50)
  mbart_target_Ns <- c(Inf, 100, 75, 50)
  warning("run_all.R: no ess calibration file found; mBART multipliers fall back to 1,0.75,0.5")
}
method_overrides <- list(
  "HierAFT.R" = list(prior_vals = 0.25),
  "mBART.R"   = list(N_target_multipliers = mbart_mult)
)
# Extra overrides for data_gen (sets number of simulated datasets to generate).
data_gen_overrides <- list(
  seed_rct    = seed_rct,  # RCT iteration seed-array generator (takes precedence over seed_method here)
  seed_rwd    = seed_rwd,  # RWD-pool seed (overrides `seed_rwd <- ...` in data_gen)
  rwd_frozen  = FALSE,     # fresh RWD per iteration (un-suffixed files); TRUE/omit = frozen pool (_fz)
  p_obs = p_obs,
  niter = 500
)
# Extra overrides for plot scripts.
plot_overrides <- list(p_obs = p_obs, mbart_target_Ns = mbart_target_Ns)
# Note: analysis files auto-infer p_obs (from data's $X columns) and niter
# (from the count of simulated data files), so neither is overridden there.

# ---- Helpers ----
format_r_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.character(x) && length(x) == 1) return(paste0('"', x, '"'))
  paste(deparse(x), collapse = " ")
}

apply_overrides <- function(text, overrides) {
  # Special-case: seed_rct (data_gen) / seed_method (analysis) rewrite the
  # `set.seed(<literal>)` seed-array generator. Per-iteration calls like
  # `set.seed(seed[ee])` are NOT touched (the regex only matches digit literals).
  for (sk in c("seed_rct", "seed_method")) {   # data_gen has seed_rct; analysis has seed_method
    if (!is.null(overrides[[sk]])) {
      text <- gsub("set\\.seed\\(\\s*\\d+\\s*\\)",
                   paste0("set.seed(", overrides[[sk]], ")"),
                   text, perl = TRUE)
      overrides[[sk]] <- NULL
      break
    }
  }

  for (name in names(overrides)) {
    val <- overrides[[name]]
    if (is.null(val)) next
    val_str <- format_r_value(val)
    pattern <- paste0("(?m)^(\\s*)\\b", name, "\\b\\s*<-\\s*[^\n]*$")
    replacement <- paste0("\\1", name, " <- ", val_str)
    text <- gsub(pattern, replacement, text, perl = TRUE)
  }
  text
}

run_methods_parallel <- function(files, overrides, n_cores = length(files)) {
  # Per-scenario parallel runner with live progress.
  # Each worker writes a tiny status file when it finishes; the main process
  # polls those files so we see methods completing in real time instead of
  # waiting silently for the slowest one. Dead workers (Stan/native crash)
  # are detected via PID liveness check so the loop can't hang.
  status_dir <- tempfile("run_status_")
  dir.create(status_dir)
  on.exit(unlink(status_dir, recursive = TRUE), add = TRUE)

  job_start <- Sys.time()
  wrapper <- function(f) {
    result <- run_one(f, overrides)
    # Write a small PLAIN-TEXT marker (status on line 1, elapsed on line 2).
    # Deliberately avoid saveRDS / the result pipe here: a forked worker that ran
    # multithreaded BART can deadlock while serializing. writeLines is low-alloc
    # (the original code used it and it worked); the parent rebuilds the result
    # row from this file below.
    writeLines(c(gsub("[\r\n]+", " ", result$status), format(result$elapsed_sec)),
               file.path(status_dir, paste0(f, ".done")))
    result
  }

  cat(sprintf("  Launching %d methods in parallel: %s\n",
              length(files), paste(files, collapse = ", ")))
  jobs <- lapply(files, function(f) parallel::mcparallel(wrapper(f), name = f))
  job_pids <- setNames(vapply(jobs, function(j) j$pid, integer(1)), files)

  pid_alive <- function(pid) {
    suppressWarnings(system2("kill", c("-0", as.character(pid)),
                             stdout = FALSE, stderr = FALSE)) == 0
  }

  pending <- files
  while (length(pending) > 0) {
    Sys.sleep(2)
    elapsed <- as.numeric(difftime(Sys.time(), job_start, units = "secs"))
    for (f in pending) {
      if (file.exists(file.path(status_dir, paste0(f, ".done")))) {
        cat(sprintf("  [%5.0fs] %-10s done\n", elapsed, f))
        pending <- setdiff(pending, f)
      } else if (!pid_alive(job_pids[[f]])) {
        cat(sprintf("  [%5.0fs] %-10s DIED (no result)\n", elapsed, f))
        pending <- setdiff(pending, f)
      }
    }
  }

  # Build results from the per-worker .done marker files (status + elapsed),
  # NOT parallel::mccollect(jobs): mccollect() does a BLOCKING read on each
  # child's result pipe and can hang forever when a forked multithreaded worker
  # (mBART / BARTv2 via OpenMP) leaves that pipe open. The polling loop above
  # already guarantees every worker finished. We deliberately do NOT kill the
  # children: force-killing a fork sibling corrupts shared connections and
  # breaks the other still-running workers ("invalid connection").
  out <- vector("list", length(files)); names(out) <- files
  for (f in files) {
    done_file <- file.path(status_dir, paste0(f, ".done"))
    out[[f]] <- if (file.exists(done_file)) {
      ln <- tryCatch(readLines(done_file, warn = FALSE), error = function(e) character(0))
      data.frame(file = f,
                 status = if (length(ln) >= 1) ln[1] else "FAILED: empty marker",
                 elapsed_sec = if (length(ln) >= 2) suppressWarnings(as.numeric(ln[2])) else NA_real_,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(file = f, status = "FAILED: worker died",
                 elapsed_sec = NA_real_, stringsAsFactors = FALSE)
    }
  }
  # Reap children without blocking on the result pipe (non-blocking).
  suppressWarnings(try(parallel::mccollect(jobs, wait = FALSE, timeout = 1),
                       silent = TRUE))
  do.call(rbind, out)
}

run_one <- function(file, overrides = list()) {
  path <- file.path(script_dir, file)
  # Layer per-file method_overrides (e.g., HierAFT's fixed prior_vals)
  # on top of the caller's overrides.
  if (exists("method_overrides", envir = globalenv()) &&
      !is.null(get("method_overrides", envir = globalenv())[[file]])) {
    overrides <- modifyList(
      overrides,
      get("method_overrides", envir = globalenv())[[file]]
    )
  }
  cat("\n----------------------------------------\n")
  cat("Running:", file, "\n")
  if (length(overrides)) {
    cat("Overrides:",
        paste(names(overrides), "=", sapply(overrides, format_r_value), collapse = ", "),
        "\n")
  }
  cat("----------------------------------------\n")
  t0 <- Sys.time()
  status <- tryCatch({
    src_text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    src_text <- apply_overrides(src_text, overrides)
    e <- new.env(parent = globalenv())
    eval(parse(text = src_text), envir = e)
    "OK"
  }, error = function(err) {
    paste("FAILED:", conditionMessage(err))
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf(">>> %s: %s (elapsed: %.1fs)\n", file, status, dt))
  data.frame(file = file, status = status, elapsed_sec = round(dt, 1),
             stringsAsFactors = FALSE)
}

# ============================================================
# Step 1: generate data
# ============================================================
if (run_data_gen) {
  cat("\n############################################################\n")
  cat("# Step 1: data generation\n")
  cat("############################################################\n")
  data_gen_file <- paste0("data_gen_p", p_obs, "_v2.R")
  for (hypo_v in hypo_values) {
    for (sc_v in sc_values) {
      for (subsc_v in subsc_values) {
        ov <- c(data_gen_overrides, common_overrides,
                list(hypo = hypo_v, sc = sc_v, subsc = subsc_v))
        run_one(data_gen_file, ov)
      }
    }
  }
}

# ============================================================
# Step 2: plot balance of generated data
# ============================================================
if (run_plot_balance) {
  cat("\n############################################################\n")
  cat("# Step 2: balance plot\n")
  cat("############################################################\n")
  run_one("plot_balance.R", plot_overrides)
}

# ============================================================
# Step 3: analysis sweep
# ============================================================
if (run_analysis) {
  cat("\n############################################################\n")
  cat("# Step 3: analysis sweep\n")
  cat("############################################################\n")
  scenarios <- list()
  for (hypo_v in hypo_values) {
    for (sc_v in sc_values) {
      for (subsc_v in subsc_values) {
        for (cor_v in cor_for_sc(sc_v)) {
          scenarios[[length(scenarios) + 1]] <- list(
            hypo  = hypo_v,
            sc    = sc_v,
            subsc = subsc_v,
            cor   = cor_v
          )
        }
      }
    }
  }

  all_results <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    s <- scenarios[[i]]
    ov <- modifyList(common_overrides, s)
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Scenario %d/%d: hypo=\"%s\", sc=%d, subsc=\"%s\", cor=%s\n",
                i, length(scenarios), s$hypo, s$sc, s$subsc, format_r_value(s$cor)))
    cat(sprintf("========================================\n"))
    if (use_parallel) {
      # Any methods listed in `serial_files` run AFTER the parallel batch.
      # Use this for methods whose native libraries don't survive a fork()
      # (e.g. RBesT/rstan-based methods on the non-surv side). Currently empty
      # for mapbart-sim-survival since all 5 methods are fork-safe.
      par_files <- setdiff(files, serial_files)
      ser_files <- intersect(files, serial_files)
      n <- if (is.null(n_cores)) length(par_files) else n_cores
      res_par <- if (length(par_files) > 0)
        run_methods_parallel(par_files, ov, n_cores = n) else NULL
      res_ser <- if (length(ser_files) > 0) {
        cat(sprintf("\n  Running %s serially (not fork-safe)...\n",
                    paste(ser_files, collapse = ", ")))
        do.call(rbind, lapply(ser_files, run_one, overrides = ov))
      } else NULL
      res <- do.call(rbind, Filter(Negate(is.null), list(res_par, res_ser)))
    } else {
      res_list <- lapply(files, run_one, overrides = ov)
      res <- do.call(rbind, res_list)
    }
    for (k in names(s)) res[[k]] <- s[[k]]
    all_results[[i]] <- res
  }
  results_all <- do.call(rbind, all_results)
  cat("\n========================================\n")
  cat("Analysis-sweep summary\n")
  cat("========================================\n")
  print(results_all)
}

# ============================================================
# Step 4: plot results
# ============================================================
if (run_plot) {
  cat("\n############################################################\n")
  cat("# Step 4: results plot\n")
  cat("############################################################\n")
  run_one("plot.R", plot_overrides)
}
