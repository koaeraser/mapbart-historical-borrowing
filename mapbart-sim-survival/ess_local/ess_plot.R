## ess_plot.R -- ESS-calibration figure (AFT-BART, survival; BART only).
##
## Loads a BART calibration .RData (ess_calibration.R, --outcome surv)
## and draws a single figure:
##
##   * Box plot of ESSbar (B values per s^2) from the AFT-BART script,
##     one box per s^2 in its S_grid.
##   * Horizontal dashed reference lines at 3.00/2.00/1.00/0.75/0.50 * N_target.
##   * Ghat(s^2) at each multiplier printed below each box.
##
## Works for both survival calibration setups (AFT-BART boxes only):
##   * ess_local survival   (ess_calibration.R; --outcome surv)
##   * mapbart-case-study-mm survival       (mapbart-case-study-mm/ess_calibration.R; files in mapbart-case-study-mm/res/)
##
## Usage:
##   Rscript ess_plot.R [--outcome gaussian|surv]
##                      [--resdir <dir>]            # where to look / write
##                      [--res <RData>]
##                      [--mult "3,2,1,0.75,0.5"]   # N_target multipliers for ref lines
##                      [--out <png>]
##
## With no overrides it auto-picks the most recent matching .RData in
## --resdir (default: this script's own res/).  For mapbart-case-study-mm, pass
## --outcome surv --resdir /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/res.  The
## output PNG defaults to the BART .RData's name with a .png extension,
## so OS / PFS / source variants never clobber each other.

## scriptDir: prefer commandArgs's --file= (the literal path Rscript was
## invoked with) so symlinked copies under sibling ess_local/ folders
## resolve scriptDir to *the symlink's folder* (e.g. mapbart-case-study-mm/ess_local/),
## not the canonical target.  normalizePath() would follow the symlink
## and send all reads/writes to the wrong res/.
scriptDir <- tryCatch({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1]))
  else dirname(sys.frames()[[1]]$ofile)
}, error = function(e) "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival/ess_local")
if (!nzchar(scriptDir) || is.na(scriptDir))
  scriptDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival/ess_local"

## Extract the calibrated (N, s^2) borrowing targets from a calibration object
## `cal` (fields $S_grid, $ESSbar, $N_target, $q).  One row per multiplier m in
## `mult_pool`, with N = round(m * N_target) and s2 = smallest s^2 with
## Ghat(m*N_target) = Pr{ESSbar <= m*N_target} >= q (the Pr-rule pick).  A
## multiplier is kept iff REACHABLE (Ghat >= q for some s^2); saturated ones
## (Ghat already >= q at the smallest s^2) are kept with s2 = smallest grid
## value.  Returns data.frame(N, s2, multiplier), sorted by N descending.
## ess_plot.R is the PRODUCER of this table; mBART.R / AFTv2.R / run_all.R read it.
ess_extract_targets <- function(cal, mult_pool) {
  empty <- data.frame(N = integer(0), s2 = numeric(0), multiplier = numeric(0))
  if (is.null(cal)) return(empty)
  S <- cal$S_grid; Nt <- cal$N_target; q <- cal$q
  rows <- lapply(sort(unique(mult_pool), decreasing = TRUE), function(m) {
    Gh <- colMeans(cal$ESSbar <= m * Nt, na.rm = TRUE)
    if (!any(Gh >= q, na.rm = TRUE)) return(NULL)   # unreachable: no valid s^2
    data.frame(N = as.integer(round(m * Nt)), s2 = S[which(Gh >= q)[1]], multiplier = m)
  })
  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(rows)) return(empty)
  rows[order(rows$N, decreasing = TRUE), , drop = FALSE]
}

## --- CLI ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  k <- grep(paste0("^", flag, "$"), args)
  if (length(k) == 0) return(default)
  args[k + 1]
}

outcome <- get_arg("--outcome", "gaussian")
stopifnot(outcome %in% c("gaussian", "surv"))

resDir <- get_arg("--resdir", file.path(scriptDir, "res"))

if (outcome == "surv") {
  ## Survival folder: the BART result files are named
  ## ess_calibration_<version>.RData (no "surv" infix any more, since
  ## the folder itself is survival-only).
  bart_pat <- "^ess_calibration_(?!gaussian_|map_)[^/]*\\.RData$"
} else {
  bart_pat <- "^ess_calibration_gaussian_.*\\.RData$"
}
## Exclude the pre-aggregate checkpoints AND the targets tables this script
## itself writes (*_tab.RData) from auto-discovery.
bart_pat <- paste0("(?!.*_(checkpoint|tab))", bart_pat)

## Optional scenario selector: restrict auto-discovery to the calibration whose
## res/ filename carries this tag (sc1E / sc2E / sc3E_cor<rho>).  Without it the
## most recent matching calibration is used (legacy behaviour).
scenario_sel <- get_arg("--scenario", NA_character_)
pick_latest_re <- function(pattern) {
  if (is.na(pattern)) return(NA_character_)
  fs <- list.files(resDir, full.names = TRUE)
  fs <- fs[grepl(pattern, basename(fs), perl = TRUE)]
  if (!is.na(scenario_sel) && nzchar(scenario_sel))
    fs <- fs[grepl(paste0("_", scenario_sel, "\\.RData$"), basename(fs))]
  if (!length(fs)) return(NA_character_)
  fs[order(file.info(fs)$mtime, decreasing = TRUE)][1]
}

bart_path <- get_arg("--res", pick_latest_re(bart_pat))

if (is.na(bart_path))
  stop("No BART .RData found in ", resDir,
       " for outcome '", outcome, "'.  Run the calibration first, or ",
       "pass --res <path> / --resdir <dir>.")

## Default output PNG: same basename as the BART .RData, .png extension,
## written alongside it.  Override with --out.
default_out <- sub("\\.RData$", ".png", basename(bart_path))
out_path    <- get_arg("--out", file.path(resDir, default_out))

cat("Outcome:    ", outcome, "\n")
cat("res dir:    ", resDir, "\n")
cat("BART .RData:", basename(bart_path), "\n")
cat("Output PNG: ", out_path, "\n")

## --- load ---------------------------------------------------------------

bart <- readRDS(bart_path)

## Scenario tag (sc1E / sc2E / sc3E_cor-0.5) parsed from the calibration's
## source RWD file (bart$canonical_file); appended to the output filename so
## plots from different scenarios never overwrite each other.
sc_tag <- if (!is.null(bart$scenario) && nzchar(bart$scenario)) bart$scenario else {
  cf <- if (!is.null(bart$canonical_file)) basename(bart$canonical_file) else ""
  m  <- regmatches(cf, regexpr("sc[0-9]+[A-Za-z]*(_cor-?[0-9.]+)?", cf))
  if (length(m)) m else ""
}
if (nzchar(sc_tag) && !grepl(sc_tag, basename(out_path), fixed = TRUE)) {
  out_path <- sub("(\\.[^.]+)$", paste0("_", sc_tag, "\\1"), out_path)
  cat("Scenario tag:", sc_tag, "-> output:", basename(out_path), "\n")
}

S_bart   <- bart$S_grid
ESSbar   <- bart$ESSbar              # B x K
Ghat     <- bart$Ghat                # length K
N_target <- bart$N_target
q_target <- bart$q
s2_star  <- bart$hat_s2_star
## endpoint (e.g. PFS/OS) and RWD sample size for the title bracket
## (guarded so older calibration files without these fields still plot)
brk_pre <- paste0(
  if (!is.null(bart$outcome)) paste0(bart$outcome, ", ") else "",
  if (!is.null(bart$n_RWD))   paste0("n=", bart$n_RWD, ", ") else "")

stopifnot(ncol(ESSbar) == length(S_bart),
          length(Ghat) == length(S_bart))

## --- plot ---------------------------------------------------------------

## Outlier handling: count points beyond the per-column Tukey whiskers
## (boxplot.stats), report counts, and clip the y-axis so they don't
## blow up the scale.  Points are not drawn (outline = FALSE), but the
## label below records how many were hidden in each column.
n_finite  <- colSums(is.finite(ESSbar))
n_out_col <- integer(ncol(ESSbar))
whisker_lo <- whisker_hi <- numeric(ncol(ESSbar))
for (j in seq_len(ncol(ESSbar))) {
  stats <- boxplot.stats(ESSbar[, j])
  n_out_col[j]  <- length(stats$out)
  whisker_lo[j] <- stats$stats[1]
  whisker_hi[j] <- stats$stats[5]
}
n_out_total <- sum(n_out_col)
n_total     <- sum(n_finite)
pct_out     <- if (n_total > 0) 100 * n_out_total / n_total else 0

cat(sprintf("Outliers (beyond Tukey whiskers) hidden: %d of %d (%.1f%%)\n",
            n_out_total, n_total, pct_out))
if (n_out_total > 0)
  cat("Per-s^2 hidden counts:",
      paste(sprintf("%s:%d", format(S_bart), n_out_col), collapse = "  "),
      "\n")

## N_target reference lines: multipliers (from --mult) + colours, defined
## up front so the y-range, annotation rows, and bottom margin all derive
## from one source.  --mult is comma/space-separated, sorted high->low;
## default 3, 2, 1, 0.75, 0.50.
mult_arg <- get_arg("--mult", "3.00,2.00,1.00,0.75,0.50")
ref_mult <- sort(unique(as.numeric(strsplit(mult_arg, "[, ]+")[[1]])),
                 decreasing = TRUE)
stopifnot(length(ref_mult) >= 1, all(is.finite(ref_mult)))
ref_h    <- ref_mult * N_target
## Colour per multiplier: a fixed palette keyed by value keeps the familiar
## scheme (3N blue, 2N green, 1N red, 0.75N orange, 0.50N goldenrod); any
## other multiplier gets a colour from a fallback palette.
.mult_pal <- c("3" = "blue", "2" = "darkgreen", "1" = "red",
               "0.75" = "orange", "0.5" = "goldenrod")
ref_col  <- unname(.mult_pal[as.character(ref_mult)])
.fallback <- c("purple", "brown", "magenta", "darkcyan", "darkorange3", "deeppink3")
.miss    <- which(is.na(ref_col))
if (length(.miss))
  ref_col[.miss] <- .fallback[((seq_along(.miss) - 1) %% length(.fallback)) + 1]
## Annotation-row positions (one Ghat row per multiplier + a removed row)
## and a bottom margin sized to fit them below the x-axis tick labels.
ghat_lines   <- 4.0 + 1.2 * (seq_along(ref_mult) - 1)
removed_line <- max(ghat_lines) + 1.4
footer_line  <- removed_line + 1.2
mar_bottom   <- footer_line + 1.3

## Combined y-range from whiskers + N_target reference lines.
y_all <- c(whisker_lo, whisker_hi, ref_h)
y_all <- y_all[is.finite(y_all)]
ylim  <- range(y_all, na.rm = TRUE)
yspan <- diff(ylim)
ylim[2] <- ylim[2] + 0.05 * yspan    # small top padding; all Ghat rows are now below

## x-positions: place every distinct s^2 at an evenly spaced index, so
## clustered s^2 values never overlap.  Tick labels still show the actual
## s^2 in full; only the spacing is uniform (not proportional to value).
all_s      <- sort(unique(round(S_bart, 8)))
pos_of     <- function(s) match(round(s, 8), all_s)
x_pos_bart <- pos_of(S_bart)
x_ticks    <- seq_along(all_s)
x_labs     <- format(all_s)
xlim       <- c(0.5, length(all_s) + 0.5)

png(out_path, width = 1400, height = 900, res = 150)
op <- par(mar = c(mar_bottom, 5.5, 3.5, 4.0))  # bottom margin sized to the
                                         # G(...N) rows + removed: k/n, all
                                         # below the rotated x-axis tick
                                         # labels; wider left margin for the
                                         # row headers

boxplot(
  ESSbar,
  at      = x_pos_bart,
  names   = format(S_bart),
  boxwex  = 0.6,
  xlim    = xlim,
  ylim    = ylim,
  xlab    = expression(s^2),
  ylab    = "ESS",
  main    = sprintf("Prior ESS calibration [%s] -- MAP-%s (boxes)\n(%sgMAP=%s, N_target=%g, B=%d, q=%g)",
                    outcome,
                    if (outcome == "surv") "AFT-BART" else "BART",
                    brk_pre, (if (!is.null(bart$map_model)) bart$map_model else bart$gmap_version), N_target, bart$B, q_target),
  outline = FALSE,                  # outliers hidden (annotated below)
  col     = "grey85",
  border  = "grey30",
  axes    = FALSE
)
axis(1, at = x_ticks, labels = x_labs, las = 2, cex.axis = 0.8)
axis(2)
box()

abline(h = ref_h, lty = 2, col = ref_col, lwd = 1.5)
## Right-edge labels: full "m x N (=value)" form, colour-matched to
## each reference line.  Replaces the corresponding legend entries.
text(x      = par("usr")[2],
     y      = ref_h,
     labels = sprintf(" %.2f x N (=%g)", ref_mult, ref_h),
     pos    = 2, cex = 0.75, col = ref_col, font = 2)

## Per-multiplier Ghat and s2_star pick (Pr-rule).
Ghat_mult <- matrix(sapply(ref_mult, function(m)
                            colMeans(ESSbar <= m * N_target, na.rm = TRUE)),
                     nrow = length(S_bart))  # K x M (force matrix shape; sapply/vapply collapse to a vector when K=1)
pick_star <- function(g) { k <- which(g >= q_target)
                           if (length(k)) S_bart[min(k)] else NA }
s2_star_mult <- vapply(seq_along(ref_mult),
                        function(mi) pick_star(Ghat_mult[, mi]), numeric(1))

## Vertical lines at each multiplier's s2_star, colour-matched.
## (Plotted before the boxes/labels so they sit under everything else.)
for (mi in seq_along(ref_mult))
  if (!is.na(s2_star_mult[mi]))
    abline(v = pos_of(s2_star_mult[mi]), lty = 3, col = ref_col[mi], lwd = 1.2)

## All annotations sit BELOW the rotated x-axis tick labels (which occupy
## lines ~1-2.5).  One colour-matched row per N_target multiplier (header
## just left of the first box), then a final "removed: k/n" row.  Rows
## and headers are generated from ref_mult so they always stay in sync
## with the reference lines above.
header_x <- xlim[1] - 0.015 * diff(xlim)
## Mark the selected s2_star cell in each Ghat row with a trailing "*".
mark_pick <- function(vals, s_pick) {
  if (is.na(s_pick)) return(vals)
  j <- which(abs(S_bart - s_pick) < 1e-8)
  if (length(j)) vals[j] <- paste0(vals[j], "*")
  vals
}
bottom_rows <- lapply(seq_along(ref_mult), function(mi)
  list(line   = ghat_lines[mi],
       vals   = mark_pick(sprintf("%.2f", Ghat_mult[, mi]), s2_star_mult[mi]),
       header = sprintf("G(%.2fN):", ref_mult[mi]),
       col    = ref_col[mi], cex = 0.75, bold = TRUE))
bottom_rows  <- c(bottom_rows, list(
  list(line   = removed_line,
       vals   = sprintf("%d/%d", n_out_col, n_finite),
       header = "removed:",
       col    = ifelse(n_out_col > 0, "grey25", "grey55"),
       cex    = 0.7, bold = FALSE)))
for (r in bottom_rows) {
  mtext(r$vals, side = 1, at = x_pos_bart, line = r$line,
        cex = r$cex, col = r$col)
  mtext(r$header, side = 1, at = header_x, line = r$line,
        cex = r$cex,
        col = if (length(unique(r$col)) == 1) r$col[1] else "grey30",
        adj = 1, font = if (r$bold) 2 else 1)
}

## Brief footer note explaining the "*" convention.
mtext("(* in each row marks the selected s2_star: smallest s^2 with Ghat >= q)",
      side = 1, line = footer_line, adj = 0, cex = 0.7, col = "grey30")

par(op)
invisible(dev.off())

cat("Wrote", out_path, "\n")

## Also print Ghat(s^2) for each N_target multiplier to the console.
ghat_df <- data.frame(s2 = S_bart,
                       medH_meanB_ESS = round(colMeans(ESSbar, na.rm = TRUE), 1))
for (mi in seq_along(ref_mult))
  ghat_df[[sprintf("Ghat_%.2fN", ref_mult[mi])]] <- round(Ghat_mult[, mi], 3)
cat("\nGhat(s^2) from BART script (one column per N_target multiplier):\n")
print(ghat_df)

## Selected s^2 for each N_target multiplier.
##   BART (Pr-rule): smallest s^2 in S_bart with Ghat(m*N_target) >= q
pick <- function(s_grid, ok) {
  k <- which(ok)
  if (length(k)) s_grid[min(k)] else NA_real_
}
sel_bart <- sapply(seq_along(ref_mult), function(mi)
                    pick(S_bart, Ghat_mult[, mi] >= q_target))

sel_df <- data.frame(
  multiplier      = ref_mult,
  target          = ref_mult * N_target,
  s2_star_BART    = sel_bart
)
cat(sprintf("\nSelected s^2_star for each N_target multiplier (q=%g for BART):\n",
            q_target))
print(sel_df)

## --- targets table (single source of truth for downstream borrowing) -------
## Write the calibrated (N, s^2) targets next to this calibration file, using
## the SAME multiplier pool drawn above (ref_mult): N = mult * N_target (from the
## calibration), s2 = the Pr-rule pick.  One row per REACHABLE multiplier (see
## ess_extract_targets) -- these are the (N, s^2) pairs mBART.R runs and whose N
## feeds AFTv2.R's rwd_w.  The N=Inf max-borrowing reference is NOT a calibration
## target and is added by mBART.R itself.
targets      <- ess_extract_targets(bart, ref_mult)
targets_path <- sub("\\.RData$", "_tab.RData", bart_path)
saveRDS(targets, targets_path)
cat(sprintf("\nWrote targets table (%d reachable N of %d multipliers) -> %s\n",
            nrow(targets), length(unique(ref_mult)), basename(targets_path)))
print(targets)
if (nrow(targets) < length(unique(ref_mult)))
  cat(sprintf("  (%d multiplier(s) dropped: unreachable -- Ghat(mult*N_target) never >= q here)\n",
              length(unique(ref_mult)) - nrow(targets)))
