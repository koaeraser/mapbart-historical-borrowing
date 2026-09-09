rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm"
outDir  <- file.path(mainDir, "res")
library(survival)

# ============================================================
# Analysis method: raw Kaplan-Meier RMST ratio (EloKRd / UCMM).
# Runs alongside AFTv2 / HierAFT / mBART under run_all.R and writes a
# result file with the same schema (delta_hat, ci_95, settings) so
# run_all.R's extract_row() reads it uniformly.  OUTCOME (PFS default,
# or OS) selects which endpoint's RMST ratio is SAVED; both endpoints
# are still printed and the KM curve is plotted.
# ============================================================
OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
cat(sprintf("=== km_analysis.R: OUTCOME = %s ===\n", OUTCOME))

# ============================================================
# km_<OUTCOME>_<tag>.png -- single-panel Kaplan-Meier curve for the
# selected OUTCOME (PFS or OS), EloKRd (red) vs UCMM (blue), no overlay.
# ============================================================
# Select the merged cohort by FILE NAME via MERGED_FILE; default is the
# primary cohort (all regimens; E-Rd excluded). The n<N> token from the file name tags the output plot.
merged_file <- Sys.getenv("MERGED_FILE",
                          unset = file.path(mainDir, "data_cleaned/merged_elokrd_ucmm_n283.RData"))
# Allow a bare filename: resolve against data_cleaned/.
if (!file.exists(merged_file) &&
    file.exists(file.path(mainDir, "data_cleaned", basename(merged_file))))
  merged_file <- file.path(mainDir, "data_cleaned", basename(merged_file))
if (!file.exists(merged_file))
  stop("Merged file not found: ", merged_file, "  (set MERGED_FILE; run data_merge.R first)")
load(merged_file)
data_tag <- regmatches(basename(merged_file), regexpr("n[0-9]+", basename(merged_file)))
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(merged_file))
cat(sprintf("km_analysis.R: merged file = %s  (tag=%s)\n", basename(merged_file), data_tag))
# Convert survival times months -> years (consistent across analyses).
merged$pfs_years <- merged$pfs_months / 12
merged$os_years  <- merged$os_months  / 12
trt_idx  <- which(merged$trt == 1)
ctrl_idx <- which(merged$trt == 0)

col_trt  <- "#c1272d"   # EloKRd
col_ctrl <- "#0000a7"   # UCMM
xmax_pfs   <- 7    # years (EloKRD follow-up extends to ~6.75 yr)
landmark   <- 5    # years -- reference bar on the curve
risk_times <- 0:xmax_pfs

draw_km_panel <- function(title, y_field, e_field, idx_t, idx_c, xmax,
                          show_xaxis = TRUE) {
  yt <- merged[[y_field]][trt_idx[idx_t]];  et <- merged[[e_field]][trt_idx[idx_t]]
  yc <- merged[[y_field]][ctrl_idx[idx_c]]; ec <- merged[[e_field]][ctrl_idx[idx_c]]
  plot(NA, NA, xlim = c(0, xmax), ylim = c(0, 1),
       xlab = if (show_xaxis) "Years" else "", ylab = "Survival", main = title,
       xaxt = if (show_xaxis) "s" else "n",
       cex.main = 0.95, cex.lab = 0.9)
  abline(h = 0.5, col = "gray80", lty = 3)
  # 5-year landmark bar
  abline(v = landmark, col = "#888888", lty = 2, lwd = 1.3)
  text(landmark, 1.0, sprintf("%g yr", landmark), col = "#888888",
       cex = 0.7, pos = 4, offset = 0.15)

  legend_text <- character(0); legend_cols <- character(0); legend_lty <- numeric(0)
  if (length(yt) >= 2) {
    lines(survfit(Surv(yt, et) ~ 1), col = col_trt, lwd = 2,
          conf.int = FALSE, mark.time = TRUE)
    legend_text <- c(legend_text, sprintf("EloKRd (n=%d)", length(idx_t)))
    legend_cols <- c(legend_cols, col_trt); legend_lty <- c(legend_lty, 1)
  }
  if (length(yc) >= 2) {
    lines(survfit(Surv(yc, ec) ~ 1), col = col_ctrl, lwd = 2,
          conf.int = FALSE, mark.time = TRUE)
    legend_text <- c(legend_text, sprintf("UCMM   (n=%d)", length(idx_c)))
    legend_cols <- c(legend_cols, col_ctrl); legend_lty <- c(legend_lty, 1)
  }
  p_lr <- NA
  if (length(yt) >= 2 && length(yc) >= 2 && sum(c(et, ec)) > 0) {
    Yall <- c(yt, yc); Eall <- c(et, ec)
    Gall <- factor(c(rep("EloKRd", length(yt)), rep("UCMM", length(yc))),
                   levels = c("EloKRd", "UCMM"))
    lr <- tryCatch(survdiff(Surv(Yall, Eall) ~ Gall), error = function(e) NULL)
    if (!is.null(lr)) p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  }
  if (!is.na(p_lr)) {
    legend_text <- c(legend_text, sprintf("log-rank p = %.3f", p_lr))
    legend_cols <- c(legend_cols, NA); legend_lty <- c(legend_lty, NA)
  }
  legend("bottomleft", legend = legend_text,
         col = legend_cols, lwd = ifelse(is.na(legend_cols), NA, 2),
         lty = legend_lty, seg.len = 1.8, bty = "n", cex = 0.78)
}

# Number-at-risk table, drawn in its own panel and aligned to the curve's
# x-axis (same xlim and left/right margins as the curve panel above).
draw_risk_table <- function(times, y_field, e_field, idx_t, idx_c, xmax) {
  yt <- merged[[y_field]][trt_idx[idx_t]];  et <- merged[[e_field]][trt_idx[idx_t]]
  yc <- merged[[y_field]][ctrl_idx[idx_c]]; ec <- merged[[e_field]][ctrl_idx[idx_c]]
  # Three stacked sub-tables: number at risk, cumulative events, cumulative censored.
  risk   <- function(y, t)       sapply(t, function(tt) sum(y >= tt))             # at risk at t
  cumcnt <- function(y, s, t, v) sapply(t, function(tt) sum(s == v & y <= tt))    # cumulative by t
  plot(NA, NA, xlim = c(0, xmax), ylim = c(0, 11), axes = FALSE, xlab = "", ylab = "")
  cn  <- 0.74
  row <- function(y, vals, col) text(times, y, vals, col = col, cex = cn)
  arm <- function(y, lab, col)  mtext(lab, side = 2, at = y, las = 1, col = col, cex = 0.66, line = 0.3)
  hdr <- function(y, txt)       text(0, y, txt, adj = c(0, 0.5), font = 2, cex = 0.72, xpd = NA)

  hdr(10.5, "Number at risk")
  row(9.5, risk(yt, times),        col_trt);  arm(9.5, "EloKRd", col_trt)
  row(8.6, risk(yc, times),        col_ctrl); arm(8.6, "UCMM",   col_ctrl)

  hdr(6.9, "Cumulative no. of events")
  row(5.9, cumcnt(yt, et, times, 1), col_trt);  arm(5.9, "EloKRd", col_trt)
  row(5.0, cumcnt(yc, ec, times, 1), col_ctrl); arm(5.0, "UCMM",   col_ctrl)

  hdr(3.3, "Cumulative no. censored")
  row(2.3, cumcnt(yt, et, times, 0), col_trt);  arm(2.3, "EloKRd", col_trt)
  row(1.4, cumcnt(yc, ec, times, 0), col_ctrl); arm(1.4, "UCMM",   col_ctrl)

  lab <- ifelse(times %% 1 == 0, formatC(times, format = "d"),
                formatC(times, format = "f", digits = 2))
  axis(1, at = times, labels = lab, cex.axis = 0.95)
  mtext("Years", side = 1, line = 1.9, cex = 0.9)
}

# ============================================================
# RMST and median survival with 95% CIs, EloKRd vs UCMM
# ============================================================
# survRM2::rmst2 gives per-arm RMST + CI and the between-arm
# difference/ratio with CI and p-value. Median survival + CI
# come from the KM fit (NA when the curve never reaches 0.5).
library(survRM2)

# RMST horizon: tau = 5 yr to match the HierAFT / mBART RMST(5yr) ratio
# estimands (so the KM row is comparable in run_all.R's table), capped at
# the shorter arm's max follow-up so rmst2 stays defined.
TAU_TARGET <- 5   # years

summarize_endpoint <- function(label, y_field, e_field) {
  time   <- merged[[y_field]]
  status <- merged[[e_field]]
  ok     <- !is.na(time) & !is.na(status)
  time   <- time[ok]; status <- status[ok]; arm <- merged$trt[ok]  # 1=EloKRd, 0=UCMM

  cat(sprintf("\n================ %s ================\n", label))

  # --- Median survival + 95% CI per arm (years) ---
  cat("\n-- Median survival (years), KM with 95% CI --\n")
  med <- summary(survfit(Surv(time, status) ~ arm))$table
  rownames(med) <- c("UCMM (arm=0)", "EloKRd (arm=1)")
  print(round(med[, c("median", "0.95LCL", "0.95UCL")], 2))

  # --- RMST per arm + between-arm difference/ratio (tau = 5 yr, capped) ---
  tau <- min(TAU_TARGET, max(time[arm == 1]), max(time[arm == 0]))
  fit <- rmst2(time = time, status = status, arm = arm, tau = tau)
  cat(sprintf("\n-- RMST (years) up to tau = %.2f yr, with 95%% CI --\n", tau))
  cat(sprintf("  EloKRd (arm=1): %.2f  (%.2f, %.2f)\n",
              fit$RMST.arm1$rmst["Est."],
              fit$RMST.arm1$rmst["lower .95"], fit$RMST.arm1$rmst["upper .95"]))
  cat(sprintf("  UCMM   (arm=0): %.2f  (%.2f, %.2f)\n",
              fit$RMST.arm0$rmst["Est."],
              fit$RMST.arm0$rmst["lower .95"], fit$RMST.arm0$rmst["upper .95"]))
  cat("\n-- RMST difference (EloKRd - UCMM) & ratio, 95% CI --\n")
  res <- fit$unadjusted.result
  print(round(res[!grepl("^RMTL", rownames(res)), , drop = FALSE], 3))

  r  <- fit$unadjusted.result["RMST (arm=1)/(arm=0)", ]   # EloKRd / UCMM (ratio)
  rd <- fit$unadjusted.result["RMST (arm=1)-(arm=0)", ]   # EloKRd - UCMM (difference)
  invisible(list(fit = fit, tau = tau,
                 delta_hat = unname(r["Est."]),
                 ci_95 = c(unname(r["lower .95"]), unname(r["upper .95"])),
                 delta_diff = unname(rd["Est."]),
                 ci_diff_95 = c(unname(rd["lower .95"]), unname(rd["upper .95"])),
                 rmst_trt_est  = c(est = unname(fit$RMST.arm1$rmst["Est."]),
                                   lo  = unname(fit$RMST.arm1$rmst["lower .95"]),
                                   hi  = unname(fit$RMST.arm1$rmst["upper .95"])),
                 rmst_ctrl_est = c(est = unname(fit$RMST.arm0$rmst["Est."]),
                                   lo  = unname(fit$RMST.arm0$rmst["lower .95"]),
                                   hi  = unname(fit$RMST.arm0$rmst["upper .95"]))))
}

km_pfs <- summarize_endpoint("PFS", "pfs_years", "pfs_status")
km_os  <- summarize_endpoint("OS",  "os_years",  "os_status")

# ============================================================
# Save the OUTCOME endpoint's RMST(5yr) ratio for run_all.R.
# Same result schema as AFTv2/HierAFT/mBART so extract_row() reads it
# uniformly; KM borrows nothing, so target_N is NA.
# ============================================================
sel <- if (OUTCOME == "OS") km_os else km_pfs
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
km_results <- list(
  method    = "KM",
  delta_hat = sel$delta_hat,
  ci_95     = sel$ci_95,
  delta_diff = sel$delta_diff,
  ci_diff_95 = sel$ci_diff_95,
  tau_rmst  = sel$tau,
  rmst_trt_est   = sel$rmst_trt_est,
  rmst_ctrl_est  = sel$rmst_ctrl_est,
  sigma_trt_est  = c(est = NA_real_, lo = NA_real_, hi = NA_real_),  # KM: nonparametric, no residual sigma
  sigma_ctrl_est = c(est = NA_real_, lo = NA_real_, hi = NA_real_),
  rmst_fit  = sel$fit,
  settings  = list(target_N = NA_integer_, tau_rmst = sel$tau,
                   estimand = "RMST ratio (EloKRd/UCMM)")
)
km_out <- file.path(outDir, sprintf("KM_results_%s_%s.RData", OUTCOME, data_tag))
saveRDS(km_results, file = km_out)
cat(sprintf("\nKM result (RMST(%.2gyr) ratio = %.3f [%.3f, %.3f]) saved to: %s\n",
            sel$tau, sel$delta_hat, sel$ci_95[1], sel$ci_95[2], km_out))

# Plot the OUTCOME's KM curves (PFS or OS), tagged like the result file:
# km_<OUTCOME>_<tag>.png  (parallel to KM_results_<OUTCOME>_<tag>.RData).
plot_y <- if (OUTCOME == "OS") "os_years"  else "pfs_years"
plot_e <- if (OUTCOME == "OS") "os_status" else "pfs_status"
km_png <- file.path(outDir, sprintf("km_%s_%s.png", OUTCOME, data_tag))
png(km_png, width = 6, height = 7.6, units = "in", res = 300)
# Stacked layout: KM curve on top, risk/events/censored tables below (shared x-axis).
layout(matrix(c(1, 2), nrow = 2), heights = c(4, 3.2))
par(mar = c(0.6, 4, 2.5, 1), mgp = c(2.3, 0.7, 0), tcl = -0.3)
draw_km_panel(OUTCOME, plot_y, plot_e,
              seq_along(trt_idx), seq_along(ctrl_idx), xmax_pfs,
              show_xaxis = FALSE)
par(mar = c(3.4, 4, 0.4, 1), mgp = c(2.3, 0.7, 0), tcl = -0.3)
draw_risk_table(risk_times, plot_y, plot_e,
                seq_along(trt_idx), seq_along(ctrl_idx), xmax_pfs)
dev.off()

cat(sprintf("Plot saved: %s\n", km_png))
