# build_table.R: results_table.csv in the layout of Table rmst-results (one row
# per method, PFS and OS side by side) with the paper's published values, the
# reruns under both codings, and the LRC-BART rows; plus tail diagnostics.
APP <- file.path(Sys.getenv("LRC_REVISION_ROOT", "revision"), "04-application")
ref <- read.csv(file.path(APP, "results/reference_summary.csv"))
lrc <- read.csv(file.path(APP, "results/lrcbart_summary.csv"))
R <- readRDS(file.path(APP, "results/lrcbart_draws.rds"))
paper <- data.frame(method = c("KM (unadjusted)", "AFT, complete pooling", "AFT, hierarchical (s2 = 0.5)", "Standard BART",
                               "MAP-AFT-BART (N_target = 30)", "MAP-AFT-BART (N_target = 22.5)", "MAP-AFT-BART (N_target = 15)"),
                    pfs_est = c(1.02, 1.03, 1.21, 1.06, 1.08, 1.09, 1.10), pfs_lo = c(0.85, 0.71, 0.73, 0.79, 0.66, 0.65, 0.65),
                    pfs_hi = c(1.22, 1.61, 2.68, 1.51, 16.4, 21.8, 32.8),
                    os_est = c(0.90, 0.90, 1.18, 0.90, 0.97, 0.98, 1.00), os_lo = c(0.79, 0.71, 0.77, 0.76, 0.74, 0.73, 0.73),
                    os_hi = c(1.03, 1.17, 1.89, 1.08, 4.77, 5.55, 7.22))
rows <- list()
add <- function(method, source, coding, estimand, p, o, extra = list(Np = NA, ep = NA, cp = NA, No = NA, eo = NA, co = NA)) {
  rows[[length(rows) + 1]] <<- data.frame(method = method, source = source, coding = coding, estimand = estimand,
    pfs_mean = p$mean, pfs_median = p$median, pfs_lo = p$lo, pfs_hi = p$hi,
    os_mean = o$mean, os_median = o$median, os_lo = o$lo, os_hi = o$hi,
    N_target_pfs = extra$Np, ess0_pfs = extra$ep, ceiling_pfs = extra$cp,
    N_target_os = extra$No, ess0_os = extra$eo, ceiling_os = extra$co)
}
for (i in seq_len(nrow(paper))) add(paper$method[i], "paper (Table rmst-results, posterior median)", "orig", "rmst",
  list(mean = NA, median = paper$pfs_est[i], lo = paper$pfs_lo[i], hi = paper$pfs_hi[i]),
  list(mean = NA, median = paper$os_est[i], lo = paper$os_lo[i], hi = paper$os_hi[i]))
labs <- c(KM = "KM (unadjusted)", `AFT-CP` = "AFT, complete pooling", `HierAFT-0.05` = "AFT, hierarchical (s2 = 0.05)",
          `HierAFT-0.5` = "AFT, hierarchical (s2 = 0.5)", BART = "Standard BART (lrcbart plain AFT-BART, 50 trees)")
for (cd in c("none", "orig", "harm")) for (m in names(labs)) for (es in c("rmst", "med")) {
  p <- ref[ref$coding == cd & ref$method == m & ref$estimand == es & ref$outcome == "PFS", ]
  o <- ref[ref$coding == cd & ref$method == m & ref$estimand == es & ref$outcome == "OS", ]
  if (nrow(p) == 0) next
  add(labs[m], "rerun (Yunxuan's scripts, or lrcbart for BART)", cd, es, p, o)
}
tgl <- c(m1 = "N_target = 30", m075 = "N_target = 22.5", m05 = "N_target = 15", f90 = "0.9 x ceiling", f50 = "0.5 x ceiling", f25 = "0.25 x ceiling")
cfgs <- unique(lrc[, c("coding", "H_f", "w")])
for (i in seq_len(nrow(cfgs))) for (tg in names(tgl)) for (es in c("rmst", "med")) {
  s <- lrc[lrc$coding == cfgs$coding[i] & lrc$H_f == cfgs$H_f[i] & lrc$w == cfgs$w[i] & lrc$target == tg & lrc$estimand == es, ]
  p <- s[s$outcome == "PFS", ]; o <- s[s$outcome == "OS", ]
  add(sprintf("LRC-BART single-arm (%s, H_f = %d, H_g = 5, w = %g)", tgl[tg], cfgs$H_f[i], cfgs$w[i]),
      "this rerun (4 chains x 2000 draws)", cfgs$coding[i], es, p, o,
      list(Np = p$N_target, ep = p$ess0, cp = p$ceiling, No = o$N_target, eo = o$ess0, co = o$ceiling))
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], function(x) signif(x, 4))
write.csv(tab, file.path(APP, "results_table.csv"), row.names = FALSE)
cat("rows:", nrow(tab), "\n")
# tail diagnostics for the reported model
cat("\nTail of the RMST ratio posterior, harmonized, H_f = 10, w = 1:\n")
for (o in c("PFS", "OS")) for (tg in c("m1", "m075", "m05")) {
  v <- as.vector(R$draws[[paste(o, "harm 10 1", tg)]]$rmst)
  cat(sprintf("%s %s: q99 %.2f q99.9 %.2f max %.2f P(>2) %.4f P(>1) %.3f\n", o, tg, quantile(v, .99), quantile(v, .999), max(v), mean(v > 2), mean(v > 1)))
}
cat("\nLRC summary columns (reported model):\n")
print(lrc[lrc$coding == "harm" & lrc$H_f == 10 & lrc$w == 1, c("outcome", "target", "estimand", "N_target", "ceiling", "s0_sq", "ess0", "ess_realized", "V_mu_f", "sigma1_sq", "tau0_sq_post", "g_sd_mean", "mean", "median", "lo", "hi", "rhat", "ess_bulk", "ess_tail")], digits = 3)
cat("\nESS map ranges (m1):\n"); for (o in c("PFS","OS")) { m <- R$maps[[paste(o, "harm 10 1 m1")]]; cat(o, round(range(m),1), " median", round(median(m),1), "\n") }
cat("\ng summary (m1): mean |g| range, sd range\n"); for (o in c("PFS","OS")) { g <- R$gsum[[paste(o, "harm 10 1 m1")]]; cat(o, round(range(abs(g$mean)),3), round(range(g$sd),3), "\n") }
