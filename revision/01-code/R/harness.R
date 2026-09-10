# harness.R: parallel, resumable simulation driver.
#
# run_sim() runs R replicates of each scenario, fits every requested method
# on the same replicate data (common random numbers across methods and
# scenarios: the data seed for replicate r is seed + r in every scenario, so
# Sc4 with delta_rwd = 0 reproduces Sc1 replicate by replicate), and writes
# one RDS per (scenario, method, replicate) under out_dir. Existing files are
# skipped, so an interrupted run resumes by calling run_sim() again with the
# same arguments. Fit failures are recorded as .err files and do not stop the
# run.
#
# summarize_sim() reads the RDS files, aggregates the per-replicate metrics
# (metrics.R) and writes summary_full.csv (every column) and
# summary_table1.csv (the paper's Table 1 layout) into out_dir.
#
# time_methods() reports mean seconds per fit per method on a few replicates.
#
# Paths: all scripts source this file relative to the 01-code root, found by
# here::here() (a .here marker sits in 01-code) or by the CODE_ROOT variable
# set by the caller.

suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
})

code_root <- function() {
  if (exists("CODE_ROOT", envir = globalenv())) return(get("CODE_ROOT", envir = globalenv()))
  if (requireNamespace("here", quietly = TRUE)) return(here::here())
  getwd()
}

source_all <- function(root = code_root()) {
  for (f in c("dgp.R", "estimands.R", "comparators.R", "metrics.R")) {
    source(file.path(root, "R", f))
  }
  invisible(root)
}

rep_seed <- function(seed, r) as.integer(seed + r)
fit_seed <- function(seed, r, method) {
  as.integer((seed + r + 7919L * (sum(utf8ToInt(method)) %% 1000L)) %% .Machine$integer.max)
}

result_path <- function(out_dir, scenario, method, r) {
  file.path(out_dir, scenario, method, sprintf("rep_%04d.rds", r))
}

fit_one <- function(dat, method, outcome, r, seed, method_args = list()) {
  fun <- get_method(method)
  set.seed(fit_seed(seed, r, method))
  t0 <- Sys.time()
  args <- c(list(dat = dat, outcome = outcome), method_args[[method]])
  fit <- do.call(fun, args)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  met <- compute_metrics(fit, dat)
  met$seconds <- secs
  list(scenario = dat$constants$scenario_id, method = method, rep = r,
       seed = rep_seed(seed, r), outcome = outcome, metrics = met,
       draws = fit$draws, surface = fit$surface, map = fit$map, ess = fit$ess,
       seconds = secs)
}

#' Run the simulation.
#' @param scenarios character vector of scenario ids (see parse_scenario()).
#' @param outcome "gaussian" or "survival".
#' @param methods character vector of registered method names.
#' @param R number of replicates.
#' @param seed base seed.
#' @param workers number of parallel workers (forked).
#' @param out_dir output directory.
#' @param method_args named list of extra argument lists per method.
#' @param gen_args extra arguments to gen_data (e.g. n_mc = 0 to skip the
#'   population truth).
run_sim <- function(scenarios, outcome, methods, R, seed = 2026, workers = 18,
                    out_dir, method_args = list(), gen_args = list(),
                    verbose = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (sc in scenarios) for (m in methods) {
    dir.create(file.path(out_dir, sc, m), recursive = TRUE, showWarnings = FALSE)
  }
  log_file <- file.path(out_dir, "run.log")
  jobs <- expand.grid(r = seq_len(R), scenario = scenarios, stringsAsFactors = FALSE)
  # Keep only jobs with at least one missing method.
  todo <- vapply(seq_len(nrow(jobs)), function(j) {
    any(!file.exists(result_path(out_dir, jobs$scenario[j], methods, jobs$r[j])))
  }, TRUE)
  jobs <- jobs[todo, , drop = FALSE]
  if (verbose) cat(sprintf("[run_sim] %d (scenario, replicate) jobs, %d methods, %d workers\n",
                           nrow(jobs), length(methods), workers))
  if (nrow(jobs) == 0) return(invisible(out_dir))
  worker <- function(j) {
    sc <- jobs$scenario[j]; r <- jobs$r[j]
    dat <- do.call(gen_data_id, c(list(id = sc, outcome = outcome, seed = rep_seed(seed, r)), gen_args))
    dat$constants$scenario_id <- sc
    for (m in methods) {
      path <- result_path(out_dir, sc, m, r)
      if (file.exists(path)) next
      res <- tryCatch(fit_one(dat, m, outcome, r, seed, method_args),
                      error = function(e) e)
      if (inherits(res, "error")) {
        writeLines(conditionMessage(res), sub("\\.rds$", ".err", path))
        cat(sprintf("%s ERROR %s %s rep %d: %s\n", format(Sys.time()), sc, m, r,
                    conditionMessage(res)), file = log_file, append = TRUE)
      } else {
        saveRDS(res, path)
      }
    }
    cat(sprintf("%s done %s rep %d\n", format(Sys.time()), sc, r), file = log_file, append = TRUE)
    TRUE
  }
  if (workers > 1) {
    mclapply(seq_len(nrow(jobs)), worker, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(jobs)), worker)
  }
  invisible(out_dir)
}

# Read every per-replicate metrics row under out_dir.
collect_results <- function(out_dir) {
  files <- list.files(out_dir, pattern = "^rep_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("no results under ", out_dir)
  rows <- lapply(files, function(f) {
    x <- readRDS(f)
    cbind(data.frame(scenario = x$scenario, method = x$method, rep = x$rep,
                     outcome = x$outcome, stringsAsFactors = FALSE), x$metrics)
  })
  rbindlist(rows, fill = TRUE)
}

summarize_sim <- function(out_dir, write = TRUE) {
  dt <- collect_results(out_dir)
  outcome <- unique(dt$outcome)
  stopifnot(length(outcome) == 1)
  agg <- aggregate_metrics(dt, outcome)
  t1 <- table1_view(agg, outcome)
  if (write) {
    fwrite(agg, file.path(out_dir, "summary_full.csv"))
    fwrite(t1, file.path(out_dir, "summary_table1.csv"))
    fwrite(dt, file.path(out_dir, "per_replicate.csv"))
  }
  list(table1 = t1, full = agg, per_replicate = dt)
}

# Mean seconds per fit per method (sequential, n_rep replicate datasets).
time_methods <- function(scenario, outcome, methods, n_rep = 3, seed = 1,
                         method_args = list()) {
  out <- lapply(methods, function(m) {
    secs <- vapply(seq_len(n_rep), function(r) {
      dat <- gen_data_id(scenario, outcome, seed = seed + r, n_mc = 0)
      dat$constants$scenario_id <- scenario
      fit_one(dat, m, outcome, r, seed, method_args)$seconds
    }, 0)
    data.frame(method = m, outcome = outcome, scenario = scenario,
               mean_seconds = mean(secs), sd_seconds = sd(secs), n_rep = n_rep)
  })
  do.call(rbind, out)
}
