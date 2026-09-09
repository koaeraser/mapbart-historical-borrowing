# Master script: full pipeline.
#
#   1. Generate data    via  data_gen_p<p_obs>_v3.R
#   2. Plot balance     via  plot_balance.R
#   3. Run analyses     via  LMv1, LMv2, LMv3, HierLM, MAP, PSCL, BARTv1, BARTv2, BARTv3, mBART
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
script_dir <- file.path(mainDir, "mapbart-sim-gaussian-single-arm")

# ---- Pipeline flags ----
run_data_gen     <- TRUE
run_plot_balance <- TRUE
run_analysis     <- TRUE
run_plot         <- TRUE

# Single-arm trial toggle (mirrors the mapbart-case-study-mm real-data application, where the
# RCT has no internal control arm and the RWD is the only control). When FALSE,
# each analysis file drops the RCT internal control (D==1 & Z==0) at load time.
# NOTE: the internal-control-only / no-borrowing methods (LMv1.R, BARTv1.R,
# MAP.R) are NOT applicable when RCT_ctrl=FALSE — drop them from `files` below.
RCT_ctrl <- FALSE

# Run analysis methods in parallel within each scenario (mclapply, fork-based).
# Set FALSE to run sequentially; falls back to sequential on Windows.
parallel_methods <- TRUE
use_parallel     <- isTRUE(parallel_methods) && .Platform$OS.type == "unix"  # parallelize methods (unix only)
n_cores          <- NULL   # NULL => use length(files); otherwise set an integer

# ---- Scenario sweep grid ----
p_obs        <- 10
sc_values    <- c(1,2)
subsc_values <- c("E")
hypo_values  <- c("alternative")  # add "null" first to calibrate thresholds
cor_for_sc   <- function(sc) if (sc == 3) 0.7 else 1

# ---- Analysis files ----
# Single-arm runs only the borrowing methods that remain valid without an
# internal control arm; the internal-control-only / no-borrowing methods
# (LMv1, BARTv1, MAP) and the others are skipped.
files <- c("LMv2.R", "HierLM.R", "BARTv2.R", "mBART.R")
# Files in `serial_files` run sequentially AFTER the parallel batch (use for
# methods whose native libraries aren't fork-safe). LMv2.R and HierLM.R both
# `library(rstan)` + `stan_model()`; rstan's RcppParallel/V8 state does not
# survive a fork(), so running them under mclapply segfaults the worker (which
# then becomes a zombie and hangs the parent's liveness poll). Run them serially
# in the main process. The parallel batch keeps only the C++/OpenMP BART methods
# (BARTv2.R, mBART.R), which are fork-safe once the C++ is pre-compiled below.
serial_files <- c("LMv2.R", "HierLM.R")

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
  RCT_ctrl    = RCT_ctrl  # single-arm toggle; no-op on files lacking the var (e.g. data_gen)
)
# Per-file overrides, layered on top of `common_overrides` and the scenario
# overrides inside `run_one()`.  Use this when methods need different values
# for the same variable.
#
# MAP.R and mBART.R no longer need `prior_vals` here: each loads its own
# calibration .RData from ess_local/res/ and iterates over the
# N_target multipliers (1.00, 0.75, 0.50 by default), labelling results
# by the integer target N rather than the raw s^2 value.
# ---- MAP-BART per-scenario, per-config s^2 calibration ----
# Each MAP-BART control config (ntree:k) reads its OWN config-matched ess
# calibration whose res/ filename carries that config's tag AND the scenario tag
# (sc1E / sc2E / sc3E_cor<rho>), written by ess_local/ess_calibration.R.  From a
# candidate multiplier pool we keep each m whose Ghat(m x N_target) crosses q at
# an INTERIOR grid point (reaches q AND is not already >= q at the smallest s^2
# -> not saturated).  The per-scenario multipliers passed to mBART.R are the
# UNION over the AVAILABLE configs of the multipliers each supports; here we
# define the helpers and the union of target Ns for plot labels.
mbart_mult_pool <- c(3, 2, 1, 0.75, 0.50)

# res/ filename tag for a scenario (mirrors the suffix ess_calibration.R appends).
ess_cal_tag <- function(sc, subsc, cor) {
  if (sc == 1 || sc == 2) paste0("sc", sc, if (!is.null(subsc)) subsc else "")
  else                    paste0("sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor)
}

# Load the config- AND sc/cor-matched gaussian calibration (EXACT tags; latest
# by mtime only among files for that exact config), or NULL.  Config tag mirrors
# ess_calibration.R / mBART.R pick_prior_for: sprintf("nt%d_k%g", ntree, k).
# NOTE: real has NO size tag (unlike surv_real); filename layout is
# ess_calibration_gaussian_<map_model>_nt<ntree>_k<k>_<sc_tag>.RData.
load_ess_cal <- function(sc, subsc, cor, ntree, k) {
  cfg_tag <- sprintf("nt%d_k%g", ntree, k)
  fs <- list.files(file.path(script_dir, "ess_local", "res"),
                   pattern = "^ess_calibration_gaussian_.*\\.RData$", full.names = TRUE)
  fs <- fs[!grepl("_checkpoint\\.RData$", fs)]
  fs <- fs[grepl(paste0("_", cfg_tag, "_", ess_cal_tag(sc, subsc, cor), "\\.RData$"),
                 basename(fs))]
  if (!length(fs)) return(NULL)
  readRDS(fs[order(file.info(fs)$mtime, decreasing = TRUE)][1])
}

# Multipliers a calibration supports (interior q-crossings only).
mbart_mult_from <- function(cal) {
  if (is.null(cal)) return(numeric(0))
  av <- Filter(function(m) {
    Gh <- colMeans(cal$ESSbar <= m * cal$N_target, na.rm = TRUE)
    any(Gh >= cal$q) && Gh[1] < cal$q
  }, mbart_mult_pool)
  sort(unlist(av), decreasing = TRUE)
}

# Candidate MAP-BART control configs (ntree:k); must match mBART.R's set / the
# override below.  plot.R overlays ALL available pairs as separate curves.
mbart_control_configs <- "5:2,10:2,50:2"
.mbart_cfgs <- Filter(Negate(is.null),
  lapply(strsplit(mbart_control_configs, "[, ]+")[[1]], function(p) {
    if (!nzchar(p)) return(NULL)
    v <- strsplit(p, ":")[[1]]
    list(ntree = as.integer(v[1]), k = as.numeric(v[2]))
  }))

# Per-SCENARIO ess summaries, unioned over the configs that HAVE a calibration
# for that sc (so no available config is starved of a multiplier it supports,
# and unavailable configs contribute nothing).
ess_mults_for_scenario <- function(sc, subsc, cor) {
  ms <- unlist(lapply(.mbart_cfgs, function(cf)
    mbart_mult_from(load_ess_cal(sc, subsc, cor, cf$ntree, cf$k))))
  sort(unique(ms), decreasing = TRUE)
}
ess_target_Ns_for_scenario <- function(sc, subsc, cor) {
  unlist(lapply(.mbart_cfgs, function(cf) {
    cal <- load_ess_cal(sc, subsc, cor, cf$ntree, cf$k)
    if (is.null(cal)) integer(0) else as.integer(round(mbart_mult_from(cal) * cal$N_target))
  }))
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

# Union of target Ns across scenarios, for the results-plot labels.  Inf = the
# max-borrowing (s^2 ~ 0) reference run mBART always adds.
.all_Ns <- unlist(lapply(scenario_grid, function(s)
  ess_target_Ns_for_scenario(s$sc, s$subsc, s$cor)))
if (length(.all_Ns)) {
  mbart_target_Ns <- c(Inf, sort(unique(.all_Ns), decreasing = TRUE))
  cat(sprintf("run_all.R: MAP-BART target Ns (union over scenarios) = %s\n",
              paste(mbart_target_Ns, collapse = ", ")))
} else {
  mbart_target_Ns <- c(Inf, 100L, 75L, 50L)
  warning("run_all.R: no sc-matched ess calibration files found; mBART target Ns fall back to Inf,100,75,50")
}
lmv4_prior_vals  <- c(0.05, 0.5) # s^2 values to sweep in HierLM.R

# LMv2 rwd_w sweep: rwd_w = 1 (full borrowing) + N_target / n_rwd_nominal for each
# finite MAP-BART N_target, so LMv2 is compared at matched effective RWD sample size.
n_rwd_nominal <- 300  # approximate RWD control count (n3 from data_gen_p10_v3.R)
.finite_Ns <- as.integer(mbart_target_Ns[is.finite(mbart_target_Ns)])
lmv2_rwd_w_vals <- sort(unique(c(1, if (length(.finite_Ns)) round(.finite_Ns / n_rwd_nominal, 4))),
                         decreasing = TRUE)
cat(sprintf("run_all.R: LMv2 rwd_w values = %s\n", paste(lmv2_rwd_w_vals, collapse = ", ")))

method_overrides <- list(
  "LMv2.R"  = list(rwd_w_vals = lmv2_rwd_w_vals),
  "HierLM.R" = list(prior_vals = lmv4_prior_vals),
  "mBART.R"  = list(mbart_control_configs = mbart_control_configs)
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
plot_overrides <- list(p_obs = p_obs, RCT_ctrl = RCT_ctrl, mbart_target_Ns = mbart_target_Ns,
                       lmv4_prior_vals = lmv4_prior_vals, lmv2_rwd_w_vals = lmv2_rwd_w_vals,
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
  # Layer in method-specific overrides (defined in the global env of
  # run_all.R) on top of whatever the caller passed.
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

  # MAP must run OUTSIDE the parallel mclapply: RBesT loads rstan, whose
  # RcppParallel thread state does not survive a fork(). Letting MAP run as
  # a sibling fork either crashes (the "did not deliver a result" warning we
  # saw before) or, worse, deadlocks the whole batch.

  # Pre-compile the BART C++ ONCE (parallel path only).  BARTv2.R and mBART.R
  # each sourceCpp() their C++; under the forked parallel runner two concurrent
  # compiles collide ("make: getcwd ..." -> "No rule to make target cmbart.o").
  # Warm the per-method sourceCpp caches here in the main process so the forked
  # workers reuse the cached .so and never invoke the compiler.  Each entry
  # mirrors a (cpp, cacheDir) pair from the method scripts exactly.
  if (use_parallel) {
    cat("Pre-compiling BART C++ (warming sourceCpp caches for the parallel batch)...\n")
    suppressPackageStartupMessages({ library(Rcpp); library(RcppEigen) })
    .precompile <- list(
      c("/wBART/cwbart.cpp", "scpp_BARTv2"),  # BARTv2.R
      c(if (isTRUE(RCT_ctrl)) "/mBART/cmbart.cpp" else "/mBART.real/cmbart.cpp",
        "scpp_mBART"),                        # mBART.R (RCT_ctrl-dependent build)
      c("/wBART/cwbart.cpp", "scpp_mBART")    # mBART.R
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
    # Per-scenario MAP-BART multipliers from the config- AND sc/cor-matched
    # calibration (union over the available configs), injected into mBART.R
    # (no-op on method files without N_target_multipliers).
    .mm <- ess_mults_for_scenario(s$sc, s$subsc, s$cor)
    if (length(.mm)) ov <- modifyList(ov, list(N_target_multipliers = .mm))
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Scenario %d/%d: hypo=\"%s\", sc=%d, subsc=\"%s\", cor=%s\n",
                i, length(scenarios), s$hypo, s$sc, s$subsc, format_r_value(s$cor)))
    cat(sprintf("========================================\n"))
    if (use_parallel) {
      # Any methods listed in `serial_files` run AFTER the parallel batch.
      # Use this for methods whose native libraries don't survive a fork()
      # (e.g. RBesT/rstan-based MAP here).
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
