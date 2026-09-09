## ess_plot_simple.R -- minimal, manuscript-ready ESS-calibration boxplot.
##
## Same data/loader as ess_plot.R, but draws ONLY the boxplot of ESSbar
## (one box per s^2): no below-axis annotation rows, no title. Outcome
## (gaussian vs surv) is auto-detected from the res/ folder; override with
## --outcome. Reference lines at the N_target multipliers are OFF by default
## (pure boxplot); turn them on with --reflines TRUE.
##
## Usage:
##   Rscript ess_plot_simple.R [--outcome gaussian|surv] [--resdir <dir>]
##                             [--res <RData>] [--mult "3,2,1,0.75,0.5"]
##                             [--reflines TRUE|FALSE] [--out <file.pdf|.png>]

## --- script dir (mirror ess_plot.R so symlinked copies resolve locally) ---
scriptDir <- tryCatch({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1]))
  else dirname(sys.frames()[[1]]$ofile)
}, error = function(e) getwd())
if (!nzchar(scriptDir) || is.na(scriptDir)) scriptDir <- getwd()

## --- CLI ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  k <- grep(paste0("^", flag, "$"), args); if (!length(k)) return(default); args[k + 1]
}
resDir <- get_arg("--resdir", file.path(scriptDir, "res"))

## Auto-detect outcome from the calibration files present (gaussian files
## carry a "_gaussian_" infix; survival files do not).
outcome <- get_arg("--outcome", NA_character_)
if (is.na(outcome)) {
  fs <- list.files(resDir)
  outcome <- if (any(grepl("^ess_calibration_gaussian_", fs))) "gaussian" else "surv"
}
stopifnot(outcome %in% c("gaussian", "surv"))

bart_pat <- if (outcome == "surv")
  "^ess_calibration_(?!gaussian_|map_)[^/]*\\.RData$" else
  "^ess_calibration_gaussian_.*\\.RData$"
bart_pat <- paste0("(?!.*_checkpoint)", bart_pat)   # skip pre-aggregate checkpoints

## Optional scenario selector (sc1E / sc2E / sc3E_cor<rho>); without it the most
## recent matching calibration is used.
scenario_sel <- get_arg("--scenario", NA_character_)
pick_latest_re <- function(pattern) {
  fs <- list.files(resDir, full.names = TRUE)
  fs <- fs[grepl(pattern, basename(fs), perl = TRUE)]
  if (!is.na(scenario_sel) && nzchar(scenario_sel))
    fs <- fs[grepl(paste0("_", scenario_sel, "\\.RData$"), basename(fs))]
  if (!length(fs)) return(NA_character_)
  fs[order(file.info(fs)$mtime, decreasing = TRUE)][1]
}
bart_path <- get_arg("--res", pick_latest_re(bart_pat))
if (is.na(bart_path))
  stop("No BART .RData found in ", resDir, " for outcome '", outcome,
       "'. Run the calibration first, or pass --res / --resdir.")

## Output: vector PDF by default (best for manuscripts); .png also supported
## via --out. Defaults to <RData basename>_simple.pdf alongside the data.
out_path <- get_arg("--out", file.path(resDir, sub("\\.RData$", "_simple.pdf", basename(bart_path))))
reflines <- toupper(get_arg("--reflines", "FALSE")) %in% c("TRUE", "T", "1", "YES")

cat("Outcome:", outcome, "| BART .RData:", basename(bart_path), "| Output:", basename(out_path), "\n")

## --- load ---------------------------------------------------------------
bart     <- readRDS(bart_path)

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
ESSbar   <- bart$ESSbar              # B x K (draws x s^2)
N_target <- bart$N_target
q_target <- bart$q
stopifnot(ncol(ESSbar) == length(S_bart))

## y-range from Tukey whiskers (outliers hidden, not drawn) + reference lines
wlo <- whi <- numeric(ncol(ESSbar))
for (j in seq_len(ncol(ESSbar))) {
  st <- boxplot.stats(ESSbar[, j]); wlo[j] <- st$stats[1]; whi[j] <- st$stats[5]
}
mult_arg <- get_arg("--mult", "3.00,2.00,1.00,0.75,0.50")
ref_mult <- sort(unique(as.numeric(strsplit(mult_arg, "[, ]+")[[1]])), decreasing = TRUE)
ref_h    <- ref_mult * N_target

## Colourblind-safe (Okabe-Ito) palette for the reference lines.
.mult_pal <- c("3" = "#0072B2", "2" = "#009E73", "1" = "#D55E00",
               "0.75" = "#E69F00", "0.5" = "#CC79A7")
ref_col  <- unname(.mult_pal[as.character(ref_mult)])
.fallback <- c("#56B4E9", "#F0E442", "#999999", "#000000")
.miss <- which(is.na(ref_col))
if (length(.miss)) ref_col[.miss] <- .fallback[((seq_along(.miss) - 1) %% length(.fallback)) + 1]

y_all <- c(wlo, whi, if (reflines) ref_h); y_all <- y_all[is.finite(y_all)]
ylim  <- range(y_all); ylim[2] <- ylim[2] + 0.04 * diff(ylim)
xlim  <- range(S_bart)
x_ticks <- sort(unique(round(S_bart, 8)))

## --- plot ---------------------------------------------------------------
if (grepl("\\.png$", out_path, ignore.case = TRUE)) {
  png(out_path, width = 1800, height = 1200, res = 300)   # 6 x 4 in @ 300 dpi
} else {
  pdf(out_path, width = 6.5, height = 4.2)                 # vector, manuscript size
}
op <- par(mar = c(4.2, 4.4, 1.0, if (reflines) 6.0 else 1.2),  # no big bottom margin, no title margin
          mgp = c(2.5, 0.7, 0), tcl = -0.3, cex.lab = 1.05, cex.axis = 0.9,
          font.lab = 1, family = "Helvetica")

boxplot(
  ESSbar, at = S_bart, names = format(S_bart),
  boxwex  = 0.6 * min(diff(sort(unique(S_bart)))),
  xlim = xlim, ylim = ylim,
  xlab = expression(italic(s)^2), ylab = "Effective sample size",
  outline = FALSE,
  col = "#CAD8E6", border = "#34495E", medcol = "#1F2D3A", medlwd = 2, whisklty = 1, staplewex = 0.5,
  axes = FALSE)
axis(1, at = x_ticks, labels = format(x_ticks), las = 2, lwd = 0, lwd.ticks = 1)
axis(2, las = 1, lwd = 0, lwd.ticks = 1)
## L-shaped axes (cleaner than a full box) for manuscript style.
box(bty = "l")

if (reflines) {
  abline(h = ref_h, lty = 2, col = ref_col, lwd = 1.3)
  ## Compact right-margin labels, colour-matched (no clutter inside the panel).
  par(xpd = NA)
  text(par("usr")[2], ref_h, labels = sprintf(" %.2gN", ref_mult),
       pos = 4, cex = 0.78, col = ref_col, font = 2)
  par(xpd = FALSE)
}

par(op); invisible(dev.off())
cat("Wrote", out_path, "\n")
