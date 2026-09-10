# make_fig_sc4.R: the borrowing map and the estimated discrepancy against X5
# under a regional shift (Sc4, R = {X5 > 2}) and a global shift (Sc5).
# The Phase 2 result files keep only the metric rows, so this script refits
# LRC-BART-100 on the first N_REP replicates of each scenario (same seeds as
# the full study: data seed 2026 + r, fit seed fit_seed(2026, r, method)),
# collects B_c(x) = P(|g(x)| < 0.5 | data) and the posterior mean of g(x) at
# every RCT control profile, and plots binned means against X5. Output:
# 05-writing/figures/fig_sc4_map.{pdf,png}; the binned values are written to
# figures/fig_sc4_map_data.csv for the consistency audit.
Sys.setenv(OMP_NUM_THREADS = "1")
REV <- Sys.getenv("LRC_REVISION_ROOT", "revision")
CODE_ROOT <- file.path(REV, "01-code")
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "fit_lrcbart.R"))
FIG <- file.path(REV, "05-writing", "figures"); dir.create(FIG, showWarnings = FALSE)
N_REP <- 40; SEED <- 2026; METHOD <- "LRC-BART-100"
scen <- c("Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d2")
delta <- c(0.5, 1, 2, 2)
fun <- get_method(METHOD)

one <- function(sc, r) {
  dat <- gen_data_id(sc, "gaussian", seed = rep_seed(SEED, r), n_mc = 0)
  set.seed(fit_seed(SEED, r, METHOD))
  fit <- fun(dat = dat, outcome = "gaussian")
  data.frame(scenario = sc, rep = r, x5 = dat$rct_ctrl$X5, map = fit$map,
             g = fit$extra$g_mean, g_sd = fit$extra$g_sd)
}
prof <- file.path(FIG, "fig_sc4_map_profiles.rds")
if (file.exists(prof)) res <- readRDS(prof) else {
  t0 <- Sys.time()
  res <- do.call(rbind, parallel::mclapply(seq_len(N_REP), function(r) {
    do.call(rbind, lapply(scen, one, r = r))
  }, mc.cores = 8))
  cat("fits done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
  saveRDS(res, prof)
}

# bins of width 0.5 on X5 over the central range (X5 ~ N(2, 1); the boundary is at 2)
brk <- seq(-0.5, 4.5, by = 0.5)
res$bin <- cut(res$x5, brk, include.lowest = TRUE)
agg <- aggregate(cbind(map, g) ~ scenario + bin, res, mean)
agg$n <- aggregate(map ~ scenario + bin, res, length)$map
agg$mid <- ((brk[-length(brk)] + brk[-1]) / 2)[as.integer(agg$bin)]
write.csv(agg, file.path(FIG, "fig_sc4_map_data.csv"), row.names = FALSE)

cols <- c("Sc4_X5_d0.5" = "#8c8c8c", "Sc4_X5_d1" = "#4477aa", "Sc4_X5_d2" = "#cc3311", "Sc5_d2" = "#228833")
lab <- c("Sc4_X5_d0.5" = expression(paste("Sc4, ", delta[2] == 0.5)),
         "Sc4_X5_d1" = expression(paste("Sc4, ", delta[2] == 1)),
         "Sc4_X5_d2" = expression(paste("Sc4, ", delta[2] == 2)),
         "Sc5_d2" = expression(paste("Sc5, ", delta[2] == 2)))
draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 1.5, 0.8), mgp = c(2.6, 0.8, 0), cex = 0.9)
  # (a) borrowing map
  plot(NA, xlim = range(brk), ylim = c(0, 1), xlab = expression(X[5]),
       ylab = expression(paste("mean ", B[c](x), ",  ", c == 0.5)), main = "(a) Borrowing map")
  abline(v = 2, lty = 3, col = "grey50")
  for (sc in scen) { a <- agg[agg$scenario == sc, ]; lines(a$mid, a$map, col = cols[sc], lwd = 2); points(a$mid, a$map, col = cols[sc], pch = 16, cex = 0.8) }
  text(2.05, 0.02, expression(paste(R == group("{", X[5] > 2, "}"))), adj = c(0, 0), cex = 0.85, col = "grey30")
  legend("topleft", legend = lab[scen], col = cols[scen], lwd = 2, pch = 16, bty = "n", cex = 0.85)
  # (b) estimated discrepancy
  plot(NA, xlim = range(brk), ylim = c(-2.3, 0.3), xlab = expression(X[5]),
       ylab = expression(paste("posterior mean ", g(x))), main = "(b) Estimated discrepancy")
  abline(v = 2, lty = 3, col = "grey50"); abline(h = 0, col = "grey80")
  for (i in seq_along(scen)) {
    sc <- scen[i]
    # truth: -delta outside R (X5 <= 2), 0 inside; Sc5: -delta everywhere
    if (sc == "Sc5_d2") segments(brk[1], -delta[i], brk[length(brk)], -delta[i], col = cols[sc], lty = 2)
    else { segments(brk[1], -delta[i], 2, -delta[i], col = cols[sc], lty = 2); segments(2, 0, brk[length(brk)], 0, col = cols[sc], lty = 2) }
    a <- agg[agg$scenario == sc, ]; lines(a$mid, a$g, col = cols[sc], lwd = 2); points(a$mid, a$g, col = cols[sc], pch = 16, cex = 0.8)
  }
  legend("bottomright", legend = c("posterior mean", "truth"), lty = c(1, 2), lwd = c(2, 1), bty = "n", cex = 0.85)
}
pdf(file.path(FIG, "fig_sc4_map.pdf"), width = 8.5, height = 3.8); draw(); dev.off()
png(file.path(FIG, "fig_sc4_map.png"), width = 8.5, height = 3.8, units = "in", res = 200); draw(); dev.off()
cat("figure written\n")
print(agg[order(agg$scenario, agg$mid), ])
