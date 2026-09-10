# Phase 2 full study, block C: the paper's single-arm design (01-code/R/
# single_arm.R). RCT treatment arm only (n_trt = 200 Gaussian; 30 and 200
# survival), RWD n = 300 the sole control source, Sc1 and Sc2, R = 500, seed
# 2026, 16 workers. Methods: LM-CP/AFT-CP, HierLM/HierAFT, BART-CP; LRC-BART
# single-arm with w = 1 (reported) and w = 0.9 (sensitivity) at N_target f90
# and 100, and LRC-BART-full (w = 1, s_0^2 at the grid minimum, the paper's
# near-complete pooling). One results directory per design; summary_table.csv
# stacks the three with a design column. Usage: Rscript run_single_arm.R
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
FULL_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
OUT_DIR <- file.path(FULL_DIR, "single_arm")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
source(file.path(CODE_ROOT, "R", "single_arm.R"))
source(file.path(FULL_DIR, "make_summary_table.R"))
register_single_arm_methods()
designs <- list(list(name = "gaussian_n200", outcome = "gaussian", n_trt = 200),
                list(name = "survival_n30", outcome = "survival", n_trt = 30),
                list(name = "survival_n200", outcome = "survival", n_trt = 200))
scenarios <- c("Sc1", "Sc2")
R <- 500; WORKERS <- 16
t0 <- Sys.time()
cat("start", format(t0), "\n")
tabs <- list(); lines <- sprintf("start %s", format(t0)); n_err <- 0
for (d in designs) {
  td <- Sys.time()
  rd <- file.path(OUT_DIR, d$name, "results")
  run_sim(scenarios, d$outcome, single_arm_methods(d$outcome), R = R, seed = 2026, workers = WORKERS,
          out_dir = rd, gen_args = single_arm_gen_args(d$n_trt))
  res <- summarize_sim(rd)
  for (f in c("summary_full.csv", "summary_table1.csv", "per_replicate.csv"))
    file.copy(file.path(rd, f), file.path(OUT_DIR, d$name, f), overwrite = TRUE)
  tab <- make_summary_table(res$full, d$outcome,
                            extra_cols = data.frame(design = d$name, outcome = d$outcome, n_trt = d$n_trt))
  tabs[[d$name]] <- tab
  errs <- list.files(rd, pattern = "\\.err$", recursive = TRUE)
  n_err <- n_err + length(errs)
  lines <- c(lines, sprintf("%s done %s wall_minutes %.1f fit_errors %d", d$name, format(Sys.time()),
                            as.numeric(Sys.time() - td, units = "mins"), length(errs)))
  cat(tail(lines, 1), "\n")
}
tab <- data.table::rbindlist(tabs, fill = TRUE)
data.table::fwrite(tab, file.path(OUT_DIR, "summary_table.csv"))
t2 <- Sys.time()
writeLines(c(lines, sprintf("summary_done %s", format(t2)),
             sprintf("wall_minutes_total %.1f", as.numeric(t2 - t0, units = "mins")),
             sprintf("R %d workers %d", R, WORKERS), sprintf("fit_errors %d", n_err)),
           file.path(OUT_DIR, "wall_time.txt"))
cat("all done", format(t2), " total wall minutes:", round(as.numeric(t2 - t0, units = "mins"), 1),
    " fit errors:", n_err, "\n")
