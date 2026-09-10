# Phase 2 full study, block A: Gaussian outcome, R = 500, seed 2026 (the dev
# seed, so replicates 1..200 are paired with dev_gaussian_v2), 16 workers.
# 12 scenarios; comparators plus LRC-BART-100/75/50, f90/f50/f25, Hg10,
# C-BART, w0. Resumable: existing RDS files are skipped.
# Usage: Rscript run_gaussian.R
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
FULL_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
OUT_DIR <- file.path(FULL_DIR, "gaussian")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
source(file.path(FULL_DIR, "make_summary_table.R"))
scenarios <- c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5",
               "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc4_X7_d1", "Sc4_X7_d2",
               "Sc5_d1", "Sc5_d2")
methods <- c(default_methods("gaussian"), "LRC-BART-100", "LRC-BART-75", "LRC-BART-50",
             "LRC-BART-f90", "LRC-BART-f50", "LRC-BART-f25", "LRC-BART-Hg10",
             "C-BART", "LRC-BART-w0")
R <- 500; WORKERS <- 16
t0 <- Sys.time()
cat("start", format(t0), "\n")
run_sim(scenarios, "gaussian", methods, R = R, seed = 2026, workers = WORKERS,
        out_dir = file.path(OUT_DIR, "results"))
t1 <- Sys.time()
cat("fits done", format(t1), " wall minutes:", round(as.numeric(t1 - t0, units = "mins"), 1), "\n")
res <- summarize_sim(file.path(OUT_DIR, "results"))
for (f in c("summary_full.csv", "summary_table1.csv", "per_replicate.csv"))
  file.copy(file.path(OUT_DIR, "results", f), file.path(OUT_DIR, f), overwrite = TRUE)
tab <- make_summary_table(res$full, "gaussian")
data.table::fwrite(tab, file.path(OUT_DIR, "summary_table.csv"))
errs <- list.files(file.path(OUT_DIR, "results"), pattern = "\\.err$", recursive = TRUE)
t2 <- Sys.time()
writeLines(c(sprintf("start %s", format(t0)), sprintf("fits_done %s", format(t1)), sprintf("summary_done %s", format(t2)),
             sprintf("wall_minutes_fits %.1f", as.numeric(t1 - t0, units = "mins")),
             sprintf("wall_minutes_total %.1f", as.numeric(t2 - t0, units = "mins")),
             sprintf("R %d workers %d scenarios %d methods %d", R, WORKERS, length(scenarios), length(methods)),
             sprintf("fit_errors %d", length(errs))),
           file.path(OUT_DIR, "wall_time.txt"))
cat("all done", format(t2), " total wall minutes:", round(as.numeric(t2 - t0, units = "mins"), 1),
    " fit errors:", length(errs), "\n")
