# analyze_full.R: post-run analysis for RESULTS_SUMMARY.md. Prints (to stdout,
# markdown) the headline tables per block, the Sc4 map and surface results,
# the fractional ladder, the go/no-go criteria at R = 500, the single-arm
# tables against the paper, and every number that changed sign or crossed a
# threshold relative to the R = 200 development run.
# Usage: Rscript analyze_full.R [gaussian_dir] (default: this directory)
suppressPackageStartupMessages(library(data.table))
FULL <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
DEV <- file.path(FULL, "..", "dev_gaussian_v2")
f3 <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
md <- function(dt, cols, digits = 3) {
  dt <- as.data.frame(dt)
  cat("| ", paste(cols, collapse = " | "), " |\n|", paste(rep("---", length(cols)), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(dt))) {
    v <- sapply(cols, function(cn) { x <- dt[i, cn]; if (is.numeric(x)) f3(x, if (cn %in% c("cov", "pow", "pow90", "map_R", "map_Rc", "map_all", "ess_capped")) 2 else digits) else as.character(x) })
    cat("| ", paste(v, collapse = " | "), " |\n", sep = "")
  }
  cat("\n")
}
headline <- c("LRC-BART-100", "LRC-BART-f90", "BART-PP", "BART-NP", "BART-CP", "C-BART")

# ---------------- Gaussian ----------------
gdir <- file.path(FULL, "gaussian")
if (file.exists(file.path(gdir, "summary_table.csv"))) {
  g <- fread(file.path(gdir, "summary_table.csv"))
  gf <- fread(file.path(gdir, "summary_full.csv"))
  lin <- g[method %in% c("LM-NP", "LM-CP", "LM-PP", "HierLM", "PSCL", "MAP-100", "MAP-75", "MAP-50")]
  best_lin <- lin[, .SD[which.min(rmse)], by = scenario][, .(scenario, best_linear = method)]
  cat("## Gaussian headline (R = 500)\n\n")
  hl <- g[method %in% headline]
  hl <- merge(hl, best_lin, by = "scenario")
  bl <- g[best_lin, on = c("scenario", method = "best_linear")]
  bl$best_linear <- bl$method
  hl <- rbind(hl, bl, fill = TRUE)
  hl[, method := factor(method, levels = c(headline, setdiff(unique(bl$method), headline)))]
  setorder(hl, scenario, method)
  md(hl, c("scenario", "method", "bias", "sd", "rmse", "cov", "pow", "pow90", "ctrl_bias", "ctrl_rmse", "ess_prior", "ess_realized", "ess_capped"))
  cat("## Sc4/Sc5 map and surface (LRC-BART-100, BART-PP)\n\n")
  s4 <- g[grepl("Sc4|Sc5", scenario) & method %in% c("LRC-BART-100", "LRC-BART-Hg10", "BART-PP", "C-BART", "LRC-BART-w0")]
  md(s4, c("scenario", "method", "bias", "rmse", "cov", "surf_rmse", "surf_rmse_R", "surf_rmse_Rc", "map_R", "map_Rc", "mapspike_R", "mapspike_Rc", "g_mean_R", "g_mean_Rc"))
  cat("## Fractional ladder and sensitivities\n\n")
  lad <- g[method %in% c("LRC-BART-50", "LRC-BART-75", "LRC-BART-100", "LRC-BART-f25", "LRC-BART-f50", "LRC-BART-f90", "LRC-BART-Hg10")]
  lad[, method := factor(method, levels = c("LRC-BART-50", "LRC-BART-75", "LRC-BART-100", "LRC-BART-f25", "LRC-BART-f50", "LRC-BART-f90", "LRC-BART-Hg10"))]
  setorder(lad, scenario, method)
  md(lad, c("scenario", "method", "bias", "rmse", "cov", "map_R", "map_Rc", "ess_prior", "ess_realized", "ess_ceiling", "ess_capped"))
  # criteria (spec Section E, as in go_nogo_report_v2)
  crit <- function(g) {
    v <- function(sc, m, col) g[scenario == sc & method == m][[col]]
    L <- "LRC-BART-100"; P <- "BART-PP"
    rows <- list()
    add <- function(sc, name, val, pass) rows[[length(rows) + 1]] <<- data.frame(scenario = sc, criterion = name, value = val, result = ifelse(pass, "PASS", "FAIL"))
    for (sc in c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5")) {
      add(sc, "abs(bias) <= 0.05", f3(v(sc, L, "bias")), abs(v(sc, L, "bias")) <= 0.05)
      add(sc, "coverage >= 0.92", f3(v(sc, L, "cov"), 3), v(sc, L, "cov") >= 0.92)
      add(sc, "RMSE <= BART-PP + 0.02", sprintf("%s vs %s", f3(v(sc, L, "rmse")), f3(v(sc, P, "rmse"))), v(sc, L, "rmse") <= v(sc, P, "rmse") + 0.02)
    }
    add("Sc1", "RMSE <= 0.90 x BART-PP", sprintf("%s vs %s", f3(v("Sc1", L, "rmse")), f3(0.9 * v("Sc1", P, "rmse"))), v("Sc1", L, "rmse") <= 0.9 * v("Sc1", P, "rmse"))
    add("Sc1", "power >= BART-PP + 0.05", sprintf("%s vs %s", f3(v("Sc1", L, "pow")), f3(v("Sc1", P, "pow"))), v("Sc1", L, "pow") >= v("Sc1", P, "pow") + 0.05)
    sc <- "Sc4_X5_d2"
    add(sc, "abs(bias) <= 0.10", f3(v(sc, L, "bias")), abs(v(sc, L, "bias")) <= 0.1)
    add(sc, "control RMSE in R <= 0.90 x BART-PP", sprintf("%s vs %s", f3(v(sc, L, "surf_rmse_R")), f3(0.9 * v(sc, P, "surf_rmse_R"))), v(sc, L, "surf_rmse_R") <= 0.9 * v(sc, P, "surf_rmse_R"))
    add(sc, "all-spike map in R >= 0.5", f3(v(sc, L, "mapspike_R")), v(sc, L, "mapspike_R") >= 0.5)
    add(sc, "all-spike map in Rc <= 0.15", f3(v(sc, L, "mapspike_Rc")), v(sc, L, "mapspike_Rc") <= 0.15)
    add(sc, "P(|g|<0.5) in R >= 0.7 (restated)", f3(v(sc, L, "map_R")), v(sc, L, "map_R") >= 0.7)
    add(sc, "P(|g|<0.5) in Rc <= 0.15 (restated)", f3(v(sc, L, "map_Rc")), v(sc, L, "map_Rc") <= 0.15)
    add(sc, "mean g within 0.25 of 0 in R", f3(v(sc, L, "g_mean_R")), abs(v(sc, L, "g_mean_R")) <= 0.25)
    add(sc, "mean g within 0.25 of -2 in Rc", f3(v(sc, L, "g_mean_Rc")), abs(v(sc, L, "g_mean_Rc") + 2) <= 0.25)
    sc <- "Sc4_X5_d1"
    add(sc, "abs(bias) <= 0.10", f3(v(sc, L, "bias")), abs(v(sc, L, "bias")) <= 0.1)
    add(sc, "overall RMSE no worse than BART-PP", sprintf("%s vs %s", f3(v(sc, L, "rmse")), f3(v(sc, P, "rmse"))), v(sc, L, "rmse") <= v(sc, P, "rmse"))
    sc <- "Sc5_d2"
    add(sc, "abs(bias) <= 0.10", f3(v(sc, L, "bias")), abs(v(sc, L, "bias")) <= 0.1)
    add(sc, "RMSE <= 1.10 x BART-NP", sprintf("%s vs %s", f3(v(sc, L, "rmse")), f3(1.1 * v(sc, "BART-NP", "rmse"))), v(sc, L, "rmse") <= 1.1 * v(sc, "BART-NP", "rmse"))
    add(sc, "all-spike map <= 0.15", f3(v(sc, L, "mapspike_Rc")), v(sc, L, "mapspike_Rc") <= 0.15)
    add(sc, "P(|g|<0.5) <= 0.15 (restated)", f3(v(sc, L, "map_Rc")), v(sc, L, "map_Rc") <= 0.15)
    do.call(rbind, rows)
  }
  cat("## Criteria at R = 500 (LRC-BART-100)\n\n")
  cf <- crit(g); cf$n <- as.character(seq_len(nrow(cf)))
  md(cf, c("n", "scenario", "criterion", "value", "result"))
  cat(sprintf("%d of %d pass\n\n", sum(cf$result == "PASS"), nrow(cf)))
  # comparison with the dev run
  if (file.exists(file.path(DEV, "summary_full.csv"))) {
    dvf <- fread(file.path(DEV, "summary_full.csv"))
    dv <- data.table(scenario = dvf$scenario, method = dvf$method, bias = dvf$ate_bias, sd = dvf$ate_sd, rmse = dvf$ate_rmse,
                     cov = dvf$ate_cover, pow = dvf$ate_power, surf_rmse_R = dvf$surf_rmse_R, map_R = dvf$map_R, map_Rc = dvf$map_Rc,
                     mapspike_R = dvf$mapspike_R, mapspike_Rc = dvf$mapspike_Rc, g_mean_R = dvf$g_mean_R, g_mean_Rc = dvf$g_mean_Rc,
                     ctrl_bias = dvf$ctrl_bias)
    cat("## Criteria flips relative to the dev run (R = 200)\n\n")
    cd <- crit(dv)
    flip <- which(cd$result != cf$result)
    if (length(flip) == 0) cat("No criterion changed result.\n\n") else {
      fl <- data.frame(n = flip, scenario = cf$scenario[flip], criterion = cf$criterion[flip], dev = cd$value[flip], full = cf$value[flip], dev_result = cd$result[flip], full_result = cf$result[flip])
      md(fl, names(fl))
    }
    cat("## Sign changes and threshold crossings (all methods, dev R = 200 vs full R = 500)\n\n")
    m <- merge(dv, g, by = c("scenario", "method"), suffixes = c("_dev", "_full"))
    ch <- list()
    chk <- function(lab, cond, a, b) { i <- which(cond); if (length(i)) ch[[length(ch) + 1]] <<- data.frame(scenario = m$scenario[i], method = m$method[i], quantity = lab, dev = f3(a[i]), full = f3(b[i])) }
    chk("ATE bias sign", sign(m$bias_dev) != sign(m$bias_full) & (abs(m$bias_dev) > 0.005 | abs(m$bias_full) > 0.005), m$bias_dev, m$bias_full)
    chk("control bias sign", sign(m$ctrl_bias_dev) != sign(m$ctrl_bias_full) & (abs(m$ctrl_bias_dev) > 0.005 | abs(m$ctrl_bias_full) > 0.005), m$ctrl_bias_dev, m$ctrl_bias_full)
    chk("|bias| crosses 0.05", (abs(m$bias_dev) <= 0.05) != (abs(m$bias_full) <= 0.05), m$bias_dev, m$bias_full)
    chk("|bias| crosses 0.10", (abs(m$bias_dev) <= 0.10) != (abs(m$bias_full) <= 0.10), m$bias_dev, m$bias_full)
    chk("coverage crosses 0.92", (m$cov_dev >= 0.92) != (m$cov_full >= 0.92), m$cov_dev, m$cov_full)
    chk("coverage crosses 0.90", (m$cov_dev >= 0.90) != (m$cov_full >= 0.90), m$cov_dev, m$cov_full)
    pp <- m[method == "BART-PP", .(scenario, rmse_pp_dev = rmse_dev, rmse_pp_full = rmse_full)]
    m2 <- merge(m, pp, by = "scenario")
    i <- which((m2$rmse_dev <= m2$rmse_pp_dev) != (m2$rmse_full <= m2$rmse_pp_full) & m2$method != "BART-PP")
    if (length(i)) ch[[length(ch) + 1]] <- data.frame(scenario = m2$scenario[i], method = m2$method[i], quantity = "RMSE vs BART-PP order", dev = sprintf("%s vs %s", f3(m2$rmse_dev[i]), f3(m2$rmse_pp_dev[i])), full = sprintf("%s vs %s", f3(m2$rmse_full[i]), f3(m2$rmse_pp_full[i])))
    chk("map_R crosses 0.5", (m$map_R_dev >= 0.5) != (m$map_R_full >= 0.5), m$map_R_dev, m$map_R_full)
    chk("map_Rc crosses 0.15", (m$map_Rc_dev <= 0.15) != (m$map_Rc_full <= 0.15), m$map_Rc_dev, m$map_Rc_full)
    if (length(ch)) md(do.call(rbind, ch), c("scenario", "method", "quantity", "dev", "full")) else cat("None.\n\n")
    # largest movements for the headline methods
    cat("## Largest bias / RMSE movements, dev to full (headline methods)\n\n")
    mh <- m[method %in% headline, .(scenario, method, bias_dev, bias_full, d_bias = bias_full - bias_dev, rmse_dev, rmse_full, d_rmse = rmse_full - rmse_dev)]
    mh <- mh[order(-pmax(abs(d_bias), abs(d_rmse)))][1:10]
    md(mh, names(mh))
  }
}

# ---------------- Survival ----------------
sdir <- file.path(FULL, "survival")
if (file.exists(file.path(sdir, "summary_table.csv"))) {
  s <- fread(file.path(sdir, "summary_table.csv"))
  lin <- s[method %in% c("AFT-NP", "AFT-CP", "AFT-PP", "HierAFT")]
  for (est in c("med", "rmst")) {
    cat(sprintf("## Survival headline, %s ratio (R = 500)\n\n", est))
    se <- s[estimand == est]
    best_lin <- lin[estimand == est][, .SD[which.min(rmse)], by = scenario][, .(scenario, method)]
    hl <- rbind(se[method %in% headline], se[best_lin, on = c("scenario", "method")])
    hl[, method := factor(method, levels = c(headline, setdiff(unique(best_lin$method), headline)))]
    setorder(hl, scenario, method)
    md(hl, c("scenario", "method", "bias", "sd", "rmse", "cov", "pow", "pow90", "ctrl_bias", "ctrl_rmse", "ess_prior", "ess_realized", "ess_capped"))
  }
  cat("## Survival Sc4/Sc5 map and surface (log scale)\n\n")
  md(s[estimand == "med" & grepl("Sc4|Sc5", scenario) & method %in% c("LRC-BART-100", "BART-PP", "C-BART")],
     c("scenario", "method", "bias", "rmse", "surf_rmse", "surf_rmse_R", "surf_rmse_Rc", "map_R", "map_Rc", "g_mean_R", "g_mean_Rc"))
  cat("## Survival fractional ladder (median ratio)\n\n")
  lad <- s[estimand == "med" & method %in% c("LRC-BART-50", "LRC-BART-75", "LRC-BART-100", "LRC-BART-f25", "LRC-BART-f50", "LRC-BART-f90")]
  md(lad, c("scenario", "method", "bias", "rmse", "cov", "ess_prior", "ess_realized", "ess_ceiling", "ess_capped"))
}

# ---------------- Single arm ----------------
adir <- file.path(FULL, "single_arm")
if (file.exists(file.path(adir, "summary_table.csv"))) {
  a <- fread(file.path(adir, "summary_table.csv"))
  cat("## Single-arm results (R = 500)\n\n")
  setorder(a, design, estimand, scenario, method)
  md(a, c("design", "estimand", "scenario", "method", "bias", "sd", "rmse", "cov", "pow", "pow90", "trt_bias", "trt_rmse", "ctrl_bias", "ctrl_rmse", "ess_prior", "ess_realized", "ess_capped"))
}
