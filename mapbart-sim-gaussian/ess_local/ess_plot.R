## ess_plot.R -- ESS-calibration figure.
##
## Loads a BART calibration .RData (ess_calibration.R) and, when
## available, the matching vanilla-MAP .RData (ess_calibration_map.R),
## and draws a single figure:
##
##   * Box plot of ESSbar (B values per s^2) from the BART script,
##     one box per s^2 in its S_grid.
##   * Horizontal dashed reference lines at 1.00/0.75/0.50 * N_target.
##   * Ghat(s^2) at each multiplier printed above each box.
##   * Overlay of the vanilla-MAP ESS curve when present (Gaussian only;
##     survival runs are BART-only).
##
## Works for all three calibration setups:
##   * ess_local Gaussian   (ess_calibration.R + ess_calibration_map.R)
##   * ess_local survival   (ess_calibration.R --outcome surv; BART boxes only)
##   * mapbart-case-study-mm survival       (mapbart-case-study-mm/ess_calibration.R; BART boxes only,
##                           files in mapbart-case-study-mm/res/)
##
## Usage:
##   Rscript ess_plot.R [--outcome gaussian|surv]
##                      [--resdir <dir>]            # where to look / write
##                      [--res <RData>] [--map <RData>]
##                      [--mult "1,0.75,0.5"]       # N_target multipliers for ref lines
##                      [--out <png>]
##
## With no overrides it auto-picks the most recent matching .RData in
## --resdir (default: this script's own res/).  For mapbart-case-study-mm, pass
## --outcome surv --resdir /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/res.  The
## output PNG defaults to the BART .RData's name with a .png extension,
## so OS / PFS / source variants never clobber each other.

## scriptDir: prefer commandArgs's --file= (the literal path Rscript was
## invoked with) so symlinked copies under sibling ess_local/ folders
## resolve scriptDir to the symlink's folder, not the canonical target.
scriptDir <- tryCatch({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1]))
  else dirname(sys.frames()[[1]]$ofile)
}, error = function(e) "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/ess_local")
if (!nzchar(scriptDir) || is.na(scriptDir))
  scriptDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/ess_local"

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
  ## Vanilla MAP exists only for the Gaussian outcome (no MAP_surv.R
  ## analysis script), so survival plots are BART-only.
  map_pat  <- NA_character_
} else {
  bart_pat <- "^ess_calibration_gaussian_.*\\.RData$"
  map_pat  <- "^ess_calibration_map_(?!surv_).*\\.RData$"
}
## Exclude the pre-aggregate checkpoints from auto-discovery.
bart_pat <- paste0("(?!.*_checkpoint)", bart_pat)

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
map_path  <- get_arg("--map",  pick_latest_re(map_pat))

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
cat("MAP  .RData:", if (is.na(map_path)) "(none)" else basename(map_path), "\n")
cat("Output PNG: ", out_path, "\n")

## --- load ---------------------------------------------------------------

bart <- readRDS(bart_path)
map  <- if (!is.na(map_path)) readRDS(map_path) else NULL

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

stopifnot(ncol(ESSbar) == length(S_bart),
          length(Ghat) == length(S_bart))

## --mult: N_target multipliers for the reference lines + Ghat rows
## (comma/space-separated, sorted high->low).  Default: 1, 0.75, 0.50.
mult_arg <- get_arg("--mult", "1.00,0.75,0.50")
ref_mult <- sort(unique(as.numeric(strsplit(mult_arg, "[, ]+")[[1]])),
                 decreasing = TRUE)
stopifnot(length(ref_mult) >= 1, all(is.finite(ref_mult)))
ref_h    <- ref_mult * N_target
## Colour per multiplier: a fixed palette keyed by value keeps the familiar
## scheme (1N red, 0.75N orange, 0.50N goldenrod, 2N green, 3N blue); any
## other multiplier gets a colour from a fallback palette.
.mult_pal <- c("3" = "blue", "2" = "darkgreen", "1" = "red",
               "0.75" = "orange", "0.5" = "goldenrod")
ref_col  <- unname(.mult_pal[as.character(ref_mult)])
.fallback <- c("purple", "brown", "magenta", "darkcyan", "darkorange3", "deeppink3")
.miss    <- which(is.na(ref_col))
if (length(.miss))
  ref_col[.miss] <- .fallback[((seq_along(.miss) - 1) %% length(.fallback)) + 1]

## Annotation-row positions (one Ghat row per multiplier + a removed row)
## and a bottom margin sized to fit them all below the x-axis tick labels.
ghat_lines   <- 4.0 + 1.2 * (seq_along(ref_mult) - 1)
removed_line <- max(ghat_lines) + 1.4
footer_line  <- removed_line + 1.2
mar_bottom   <- footer_line + 1.3

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

## Combined y-range from whiskers + N_target reference lines + map curve.
y_all <- c(whisker_lo, whisker_hi, ref_h,
            if (!is.null(map)) map$ESS else NULL)
y_all <- y_all[is.finite(y_all)]
ylim  <- range(y_all, na.rm = TRUE)
yspan <- diff(ylim)
ylim[2] <- ylim[2] + 0.05 * yspan    # small top padding; all Ghat rows are now below

## Combined x-positions: place every distinct s^2 (from either method) at an
## evenly spaced index, so clustered s^2 values never overlap.  Tick labels
## still show the actual s^2 in full; only the spacing is uniform (not
## proportional to value).  pos_of() maps an s^2 to its slot.
x_all      <- c(S_bart, if (!is.null(map)) map$S_grid else NULL)
all_s      <- sort(unique(round(x_all, 8)))
pos_of     <- function(s) match(round(s, 8), all_s)
x_pos_bart <- pos_of(S_bart)
x_ticks    <- seq_along(all_s)
x_labs     <- format(all_s)
xlim       <- c(0.5, length(all_s) + 0.5)

png(out_path, width = 1400, height = 900, res = 150)
op <- par(mar = c(mar_bottom, 5.5, 3.5, 4.0))  # bottom margin sized to the
                                         # G(...N) rows + removed: k/n, all below
                                         # the rotated x-axis tick labels; wider
                                         # left margin for the row headers

boxplot(
  ESSbar,
  at      = x_pos_bart,
  names   = format(S_bart),
  boxwex  = 0.6,
  xlim    = xlim,
  ylim    = ylim,
  xlab    = expression(s^2),
  ylab    = "ESS",
  main    = sprintf("Prior ESS calibration [%s] -- MAP-%s (boxes) vs vanilla MAP (line)\n(gMAP=%s, N_target=%g, B=%d, q=%g)",
                    outcome,
                    if (outcome == "surv") "AFT-BART" else "BART",
                    (if (!is.null(bart$map_model)) bart$map_model else bart$gmap_version), N_target, bart$B, q_target),
  outline = FALSE,                  # outliers hidden (annotated below)
  col     = "grey85",
  border  = "grey30",
  axes    = FALSE
)
axis(1, at = x_ticks, labels = x_labs, las = 2, cex.axis = 0.8)
axis(2)
box()

## N_target reference lines (multipliers/colours set above from --mult).
abline(h = ref_h, lty = 2, col = ref_col, lwd = 1.5)
## Right-edge labels: full "m x N (=value)" form, colour-matched to
## each reference line.  Replaces the corresponding legend entries.
text(x      = par("usr")[2],
     y      = ref_h,
     labels = sprintf(" %.2f x N (=%g)", ref_mult, ref_h),
     pos    = 2, cex = 0.75, col = ref_col, font = 2)

## Per-multiplier Ghat and s2_star pick (Pr-rule).
Ghat_mult <- sapply(ref_mult, function(m)
                     colMeans(ESSbar <= m * N_target, na.rm = TRUE))  # K x M
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
## just left of the first box), then a final "removed: k/n" row.  Rows are
## generated from ref_mult so they stay in sync with the reference lines.
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
bottom_rows <- c(bottom_rows, list(
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

## Overlay vanilla-MAP ESS curve, if available.
if (!is.null(map)) {
  lines(pos_of(map$S_grid), map$ESS, type = "b", pch = 17, lwd = 2, col = "blue")
}

## Legend: per-multiplier entries are gone -- the colour-matched
## right-edge "m x N (=value)" labels identify the reference lines, and
## the "*" in each Ghat row marks the selected s2_star cell.  Only the
## vanilla-MAP overlay (if any) still needs a legend entry.
if (!is.null(map))
  legend("topright", legend = "vanilla MAP ESS",
         lty = 1, col = "blue", pch = 17, bty = "n", cex = 0.9)
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

## Selected s^2 for each method x each N_target multiplier.
##   BART (Pr-rule): smallest s^2 in S_bart with Ghat(m*N_target) >= q
##   vanilla MAP:    smallest s^2 in S_map  with ESS <= m * N_target
pick <- function(s_grid, ok) {
  k <- which(ok)
  if (length(k)) s_grid[min(k)] else NA_real_
}
sel_bart <- sapply(seq_along(ref_mult), function(mi)
                    pick(S_bart, Ghat_mult[, mi] >= q_target))
sel_map <- if (!is.null(map)) {
  sapply(ref_mult, function(m)
          pick(map$S_grid, map$ESS <= m * N_target))
} else {
  rep(NA_real_, length(ref_mult))
}

sel_df <- data.frame(
  multiplier      = ref_mult,
  target          = ref_mult * N_target,
  s2_star_BART    = sel_bart,
  s2_star_MAP     = sel_map
)
cat(sprintf("\nSelected s^2_star for each method x N_target multiplier (q=%g for BART):\n",
            q_target))
print(sel_df)
if (!is.null(map)) {
  cat("\nvanilla-MAP ESS(s^2):\n")
  print(data.frame(s2 = map$S_grid, ESS = round(map$ESS, 2)))
}
