# Phase 2 development go/no-go run v2 (after Phase 2b tuning): Gaussian outcome, R = 200, seed 2026,
# 16 workers. Scenarios and methods per method_spec.md Section E (CAHB and
# LRC-BART-Hg5 excluded from the development run). Output: one RDS per
# (scenario, method, replicate) under results/, and summarize_sim() CSVs
# (summary_full.csv, summary_table1.csv, per_replicate.csv) in this directory.
# Usage: Rscript run_dev_v2.R   (resumable: existing RDS files are skipped)
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
OUT_DIR <- file.path(dirname(CODE_ROOT), "02-validation", "dev_gaussian_v2")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
scenarios <- c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5",
               "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc4_X7_d1", "Sc4_X7_d2",
               "Sc5_d1", "Sc5_d2")
methods <- c(default_methods("gaussian"), "LRC-BART-100", "LRC-BART-75", "LRC-BART-50", "C-BART", "LRC-BART-w0", "LRC-BART-f90", "LRC-BART-f50")
t0 <- Sys.time()
cat("start", format(t0), "\n")
run_sim(scenarios, "gaussian", methods, R = 200, seed = 2026, workers = 16,
        out_dir = file.path(OUT_DIR, "results"))
t1 <- Sys.time()
cat("fits done", format(t1), " wall minutes:", round(as.numeric(t1 - t0, units = "mins"), 1), "\n")
res <- summarize_sim(file.path(OUT_DIR, "results"))
for (f in c("summary_full.csv", "summary_table1.csv", "per_replicate.csv"))
  file.copy(file.path(OUT_DIR, "results", f), file.path(OUT_DIR, f), overwrite = TRUE)
t2 <- Sys.time()
writeLines(c(sprintf("start %s", format(t0)), sprintf("fits_done %s", format(t1)), sprintf("summary_done %s", format(t2)),
             sprintf("wall_minutes_fits %.1f", as.numeric(t1 - t0, units = "mins")),
             sprintf("wall_minutes_total %.1f", as.numeric(t2 - t0, units = "mins"))),
           file.path(OUT_DIR, "wall_time.txt"))
cat("all done", format(t2), " total wall minutes:", round(as.numeric(t2 - t0, units = "mins"), 1), "\n")
