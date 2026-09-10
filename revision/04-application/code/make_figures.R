# make_figures.R: figures for the application (reported model: harmonized
# coding, H_f = 10, w = 1). (a) local ESS map ESS(x) at the 30 EloKRd profiles
# against each of the six covariates; (b) posterior of the 5-year RMST ratio
# under the three calibrated targets; (c) posterior summary of the discrepancy
# g(x) at the 30 profiles (a prior draw in single-arm mode). PDF and PNG.
APP <- file.path(Sys.getenv("LRC_REVISION_ROOT", "revision"), "04-application")
source(file.path(APP, "code/prep_data.R"))
R <- readRDS(file.path(APP, "results/lrcbart_draws.rds"))
FIG <- file.path(APP, "figures")
dev_both <- function(name, w, h, expr) {
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h)
    else png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 200)
    eval(expr); dev.off()
  }
}
cols <- c(PFS = "#c1272d", OS = "#0000a7")
dat <- load_mm("PFS"); X <- dat$X_harm[dat$trt_idx, ]
covs <- c(age = "Age (years)", male = "Male", race = "Race", hispanic = "Hispanic or Latino",
          high_risk_cyto = "High-risk cytogenetics", asct = "ASCT off protocol")
race <- ifelse(X$race_Black == 1, "Black", ifelse(X$race_Other == 1, "Other", "White"))

# (a) ESS map, both endpoints, target m = 1
dev_both("fig_ess_map", 9, 6, quote({
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
  for (v in names(covs)) {
    xv <- if (v == "race") factor(race, levels = c("White", "Black", "Other")) else X[[v]]
    ylim <- range(unlist(lapply(c("PFS", "OS"), function(o) R$maps[[paste(o, "harm 10 1 m1")]])))
    ylim <- c(0, max(ylim) * 1.05)
    if (v == "age") {
      plot(NA, xlim = range(xv), ylim = ylim, xlab = covs[v], ylab = "ESS(x)", main = covs[v])
      for (o in c("PFS", "OS")) points(xv, R$maps[[paste(o, "harm 10 1 m1")]], pch = 19, col = cols[o])
    } else {
      xf <- if (is.factor(xv)) xv else factor(xv, levels = c(0, 1), labels = c("No", "Yes"))
      plot(NA, xlim = c(0.5, nlevels(xf) + 0.5), ylim = ylim, xaxt = "n", xlab = "", ylab = "ESS(x)", main = covs[v])
      axis(1, at = seq_len(nlevels(xf)), labels = levels(xf))
      for (o in c("PFS", "OS")) points(jitter(as.numeric(xf), 0.4) + ifelse(o == "PFS", -0.12, 0.12),
                                       R$maps[[paste(o, "harm 10 1 m1")]], pch = 19, col = cols[o])
    }
    if (v == "age") legend("topright", c("PFS", "OS"), pch = 19, col = cols, bty = "n")
  }
}))

# (b) posterior of the RMST ratio under the three targets
dev_both("fig_rmst_posterior", 9, 4, quote({
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
  tg <- c(m1 = "N_target = 30", m075 = "N_target = 22.5", m05 = "N_target = 15")
  lty <- c(m1 = 1, m075 = 2, m05 = 3)
  for (o in c("PFS", "OS")) {
    dl <- lapply(names(tg), function(t) density(as.vector(R$draws[[paste(o, "harm 10 1", t)]]$rmst), from = 0.5, to = 2))
    plot(NA, xlim = c(0.5, 1.8), ylim = c(0, max(sapply(dl, function(d) max(d$y)))), xlab = "5-year RMST ratio (EloKRd / UCMM)",
         ylab = "Posterior density", main = o)
    abline(v = 1, col = "gray70", lty = 3)
    for (i in seq_along(dl)) lines(dl[[i]], col = cols[o], lwd = 2, lty = lty[i])
    legend("topright", tg, lty = lty, lwd = 2, col = cols[o], bty = "n", cex = 0.85)
  }
}))

# (c) discrepancy g(x) at the 30 profiles (prior draw in single-arm mode)
dev_both("fig_g_profiles", 9, 4, quote({
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.3, 0.7, 0))
  for (o in c("PFS", "OS")) {
    gs <- R$gsum[[paste(o, "harm 10 1 m1")]]; ord <- order(X$age)
    plot(NA, xlim = c(1, 30), ylim = range(c(gs$lo, gs$hi)) * 1.05, xlab = "EloKRd profile (ordered by age)",
         ylab = "g(x), log-time scale", main = o)
    abline(h = 0, col = "gray70", lty = 3)
    segments(1:30, gs$lo[ord], 1:30, gs$hi[ord], col = cols[o])
    points(1:30, gs$mean[ord], pch = 19, col = cols[o], cex = 0.8)
  }
}))
cat("figures written\n")
