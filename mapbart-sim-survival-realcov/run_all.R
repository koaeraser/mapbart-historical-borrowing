# Master script: full pipeline.
#
#   1. Generate data    via  data_gen_p<p_obs>_v3.R
#   2. Plot balance     via  plot_balance.R
#   3. Run analyses     via  AFTv1, AFTv2, AFTv3, HierAFT, BARTv1, BARTv2, BARTv3, mBART
#   4. Plot results     via  plot.R  (median-ratio left + RMST(tau)-ratio right, unified)
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
script_dir <- file.path(mainDir, "mapbart-sim-survival-realcov")

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

# ---- Sample-size option (single switch; tags all data/result/plot files) ----
#   "default" -> RCT trt 200, RWD 300, files tagged _n200
#   "small"   -> RCT trt 30,  RWD 200, files tagged _n30
size_option  <- "default"

# ---- Scenario sweep grid ----
p_obs        <- 10
sc_values    <- c(1,2)    
subsc_values <- c("E")
hypo_values  <- c("alternative")  # add "null" first to calibrate thresholds
cor_for_sc   <- function(sc) if (sc == 3) 0.7 else 1

# ---- Analysis files ----
files <- c("AFTv2.R", "HierAFT.R", "BARTv2.R", "mBART.R")  # single-arm borrowing methods
# Files in `serial_files` run sequentially AFTER the parallel batch (use for
# methods whose native libraries aren't fork-safe). AFTv2.R and HierAFT.R both
# `library(rstan)` + `stan_model()`; rstan's RcppParallel/V8 state does not
# survive a fork(), so running them under mclapply segfaults the worker (which
# then becomes a zombie and hangs the parent's liveness poll). Run them serially.
# The parallel batch keeps only the C++/OpenMP BART methods (BARTv2.R, mBART.R),
# which are fork-safe once the C++ is pre-compiled.
serial_files <- c("AFTv2.R", "HierAFT.R")

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
  size_option = size_option,  # single switch set above; tags files _n200/_n30 across data/methods
  rmst_control_sigma = "own"  # control-arm RMST sigma: "trt"|"own"|"adaptive" (own = control's own sigma-hat; avoids inflated trt sigma at small n)
)
# Per-file overrides, layered on top of common_overrides and the scenario
# overrides inside `run_one()`.
#
# mBART.R and AFTv2.R now SELF-SERVE their borrowing targets from the (N, s^2)
# tables that ess_local/ess_plot.R writes next to each calibration
# (ess_calibration_<version>_tab.RData).  run_all.R no longer extracts the
# Pr-rule / interior-q-crossing here -- it only READS those tables to label
# plot.R (mbart_target_Ns / aftv2_rwd_w_vals below).  Run ess_plot.R on every
# calibration BEFORE this sweep; configs without a table are skipped by mBART.R.
.size_tag  <- if (size_option == "small") "n30" else "n200"

# res/ filename tag for a scenario (mirrors the suffix ess_calibration.R appends
# AFTER the size tag).
ess_cal_tag <- function(sc, subsc, cor) {
  if (sc == 1 || sc == 2) paste0("sc", sc, if (!is.null(subsc)) subsc else "")
  else                    paste0("sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor)
}

# Candidate MAP-BART control configs (ntree:k); must match mBART.R's set / the
# override below.  plot.R overlays ALL available pairs as separate curves.
mbart_control_configs <- "5:2,10:2,50:2"

# Per-SCENARIO target Ns, read (not recomputed) from the ess_plot.R targets
# tables and unioned over the control configs available for that scenario.
tgt_Ns_for_scenario <- function(sc, subsc, cor) {
  fs <- list.files(file.path(script_dir, "ess_local", "res"), full.names = TRUE,
                   pattern = paste0("^ess_calibration_.*_", .size_tag,
                                    "_nt[0-9]+_k[0-9.]+_", ess_cal_tag(sc, subsc, cor),
                                    "_tab\\.RData$"))
  if (!length(fs)) return(integer(0))
  sort(unique(unlist(lapply(fs, function(f) as.integer(readRDS(f)$N)))), decreasing = TRUE)
}

# Scenario grid (shared by the analysis sweep and the plot's target-N labels).
build_scenarios <- function() {
  out <- list()
  for (hypo_v in hypo_values)
    for (sc_v in sc_values)
      for (subsc_v in subsc_values)
        for (cor_v in cor_for_sc(sc_v))
          out[[length(out) + 1]] <- list(hypo = hypo_v, sc = sc_v,
                                          subsc = subsc_v, cor = cor_v)
  out
}
scenario_grid <- build_scenarios()

# Union of target Ns across scenarios, for the results-plot labels ONLY (the
# analysis scripts self-serve from the tables).  Read from the ess_plot.R
# targets tables; Inf = the max-borrowing (s^2 ~ 0) reference mBART always adds.
.all_Ns <- unlist(lapply(scenario_grid, function(s)
  tgt_Ns_for_scenario(s$sc, s$subsc, s$cor)))
if (length(.all_Ns)) {
  mbart_target_Ns <- c(Inf, sort(unique(.all_Ns), decreasing = TRUE))
  cat(sprintf("run_all.R: MAP-BART target Ns (%s, union over scenarios, from targets tables) = %s\n",
              .size_tag, paste(mbart_target_Ns, collapse = ", ")))
} else {
  mbart_target_Ns <- c(Inf)
  warning(sprintf(paste0("run_all.R: no %s targets tables (ess_calibration_*_tab.RData) found -- ",
                  "plot labels will show only the N=Inf reference.  Run ess_local/ess_plot.R on the ",
                  "calibrations first."), .size_tag))
}
lmv4_prior_vals <- c(0.05, 0.5)  # s^2 values to sweep in HierAFT.R

# AFTv2 rwd_w labels for the plot: rwd_w = 1 + N / n_rwd_nominal for each finite
# MAP-BART target N (global union over scenarios).  AFTv2.R computes its own
# per-scenario rwd_w the same way; this is just the union for plot.R's file scan.
n_rwd_nominal <- 300  # approximate RWD control count (n3 from data_gen_p10_v3.R)
.finite_Ns <- as.integer(mbart_target_Ns[is.finite(mbart_target_Ns)])
aftv2_rwd_w_vals <- sort(unique(c(1, if (length(.finite_Ns)) round(.finite_Ns / n_rwd_nominal, 4))),
                          decreasing = TRUE)
cat(sprintf("run_all.R: AFTv2 rwd_w values (plot labels) = %s\n", paste(aftv2_rwd_w_vals, collapse = ", ")))

# mBART.R self-serves prior_vals/target_Ns from its config/scenario-matched
# targets table; AFTv2.R self-serves rwd_w_vals from the scenario's tables.  So
# neither is injected here -- only the control-config sweep (mBART) is overridden.
method_overrides <- list(
  "HierAFT.R" = list(prior_vals = lmv4_prior_vals),
  "mBART.R"  = list(mbart_control_configs = mbart_control_configs)
)
# Extra overrides for data_gen (sets number of simulated datasets to generate).
data_gen_overrides <- list(
  seed_rct    = seed_rct,  # RCT iteration seed-array generator (takes precedence over seed_method here)
  seed_rwd    = seed_rwd,  # RWD-pool seed (overrides `seed_rwd <- ...` in data_gen)
  rwd_frozen  = FALSE,     # fresh RWD per iteration (un-suffixed files); TRUE/omit = frozen pool (_fz)
  p_obs = p_obs,
  niter = 100
) 
# Extra overrides for plot scripts.
plot_overrides <- list(p_obs = p_obs, size_option = size_option,
                       mbart_target_Ns = mbart_target_Ns,
                       aftv4_prior_vals = lmv4_prior_vals, aftv2_rwd_w_vals = aftv2_rwd_w_vals,
                       n_rwd_nominal = n_rwd_nominal,
                       mbart_control_configs = mbart_control_configs)
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
    # A crashed forked worker becomes a <defunct> zombie until the parent reaps
    # it; `kill -0` still succeeds on a zombie, so check the process state and
    # treat Z (zombie) as dead -- otherwise the poll below hangs forever.
    st <- suppressWarnings(system2("ps", c("-o", "state=", "-p", as.character(pid)),
                                   stdout = TRUE, stderr = FALSE))
    st <- trimws(paste(st, collapse = ""))
    nzchar(st) && !startsWith(st, "Z")
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
  data_gen_file <- paste0("data_gen_p", p_obs, "_v3.R")
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

  # Pre-compile the BART C++ ONCE (parallel path only).  BARTv2.R and mBART.R
  # each sourceCpp() their C++; under the forked parallel runner two concurrent
  # compiles collide (make: getcwd ... -> "cannot open the connection").  Warm
  # the per-method sourceCpp caches here in the main process so the forked
  # workers reuse the cached .so and never invoke the compiler.  Each entry
  # mirrors a (cpp, cacheDir) pair from the method scripts exactly.
  if (use_parallel) {
    cat("Pre-compiling BART C++ (warming sourceCpp caches for the parallel batch)...\n")
    suppressPackageStartupMessages({ library(Rcpp); library(RcppEigen) })
    .precompile <- list(
      c("/aBART/cabart.cpp",      "scpp_BARTv2"),  # BARTv2.R
      c("/mBART.real/cmbart.cpp", "scpp_mBART"),   # mBART.R
      c("/aBART/cabart.cpp",      "scpp_mBART")    # mBART.R
    )
    for (.pc in .precompile) {
      .cd <- file.path(tempdir(), .pc[2])
      dir.create(.cd, showWarnings = FALSE, recursive = TRUE)
      Rcpp::sourceCpp(paste0(mainDir, .pc[1]), cacheDir = .cd, env = new.env())
    }
    cat("Pre-compile done.\n")
  }

  scenarios <- scenario_grid

  all_results <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    s <- scenarios[[i]]
    ov <- modifyList(common_overrides, s)
    # mBART.R / AFTv2.R read their per-scenario targets from the ess_plot.R
    # tables themselves (matched on the sc/subsc/cor overrides above), so no
    # N_target_multipliers / rwd_w injection is needed here.
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Scenario %d/%d: hypo=\"%s\", sc=%d, subsc=\"%s\", cor=%s\n",
                i, length(scenarios), s$hypo, s$sc, s$subsc, format_r_value(s$cor)))
    cat(sprintf("========================================\n"))
    if (use_parallel) {
      # Any methods listed in `serial_files` run AFTER the parallel batch.
      # Use this for methods whose native libraries don't survive a fork()
      # (e.g. RBesT/rstan-based methods on the non-surv side). Currently empty
      # for mapbart-sim-survival-realcov since all 5 methods are fork-safe.
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
