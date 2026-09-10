# CAHB comparator run: Gaussian outcome, Sc1, Sc4_X5_d1, Sc4_X5_d2, Sc5_d1,
# Sc5_d2, R = 200, seed 2026 (so replicate r has the same data seed 2026 + r
# as the full study in ../full/gaussian and the rows are paired with it),
# 12 workers. Methods: CAHB (primary), CAHB-p10, CAHB-lit (see
# 01-code/R/cahb.R). Resumable. Usage: OMP_NUM_THREADS=1 Rscript run_cahb.R
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "..", "..", "01-code"))
HERE <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "cahb.R"))
source(file.path(HERE, "..", "full", "make_summary_table.R"))
scenarios <- c("Sc1", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d1", "Sc5_d2")
methods <- c("CAHB", "CAHB-p10", "CAHB-lit")
R <- 200; WORKERS <- 12
t0 <- Sys.time(); cat("start", format(t0), "\n")
run_sim(scenarios, "gaussian", methods, R = R, seed = 2026, workers = WORKERS,
        out_dir = file.path(HERE, "results"))
t1 <- Sys.time()
cat("fits done", format(t1), " wall minutes:", round(as.numeric(t1 - t0, units = "mins"), 2), "\n")
res <- summarize_sim(file.path(HERE, "results"))
tab <- make_summary_table(res$full, "gaussian")
data.table::fwrite(tab, file.path(HERE, "summary_table_cahb_only.csv"))
errs <- list.files(file.path(HERE, "results"), pattern = "\\.err$", recursive = TRUE)
t2 <- Sys.time()
writeLines(c(sprintf("start %s", format(t0)), sprintf("fits_done %s", format(t1)), sprintf("summary_done %s", format(t2)),
             sprintf("wall_minutes_fits %.2f", as.numeric(t1 - t0, units = "mins")),
             sprintf("wall_minutes_total %.2f", as.numeric(t2 - t0, units = "mins")),
             sprintf("R %d workers %d scenarios %d methods %d", R, WORKERS, length(scenarios), length(methods)),
             sprintf("fit_errors %d", length(errs))),
           file.path(HERE, "wall_time.txt"))
cat("all done", format(t2), " fit errors:", length(errs), "\n")
