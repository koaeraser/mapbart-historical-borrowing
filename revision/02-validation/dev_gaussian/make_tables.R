# Builds the per-scenario tables and the go/no-go criterion table for
# go_nogo_report.md from summary_full.csv (written by run_dev.R).
# Usage: Rscript make_tables.R  -> tables.md, criteria.csv in this directory.
suppressPackageStartupMessages(library(data.table))
DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(DIR) == 0) DIR <- "."
full <- fread(file.path(DIR, "summary_full.csv"))
full[, rmse_full := ate_rmse]

lin <- c("LM-NP", "LM-CP", "LM-PP", "PSCL", "MAP-100", "MAP-75", "MAP-50", "HierLM")
bart <- c("BART-NP", "BART-CP", "BART-PP")
lrc <- c("C-BART", "LRC-BART-w0", "LRC-BART-50", "LRC-BART-75", "LRC-BART-100")
scen <- c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5",
          "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc4_X7_d1", "Sc4_X7_d2", "Sc5_d1", "Sc5_d2")
f3 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 3))
f2 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 2))
f1 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 1))

out <- character(0)
for (sc in scen) {
  d <- full[scenario == sc]
  best_lin <- d[method %in% lin][order(ate_rmse)][1, method]
  rows <- d[method %in% c(bart, best_lin, lrc)]
  rows[, ord := match(method, c(best_lin, bart, lrc))]; rows <- rows[order(ord)]
  sc4 <- grepl("^Sc4", sc)
  hdr <- c("method", "bias", "SD", "RMSE", "cov", "power")
  if (sc4) hdr <- c(hdr, "RMSE_R", "RMSE_Rc", "map_R", "map_Rc", "spike_R", "spike_Rc", "g_R", "g_Rc", "ESS0", "ESSreal", "ceiling", "capped")
  else hdr <- c(hdr, "map", "g", "ESS0", "ESSreal", "ceiling", "capped")
  out <- c(out, sprintf("### %s (n_rep = %d; best linear method: %s)", sc, rows$n_rep[1], best_lin), "",
           paste0("| ", paste(hdr, collapse = " | "), " |"),
           paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"))
  for (i in seq_len(nrow(rows))) {
    r <- rows[i]
    cells <- c(r$method, f3(r$ate_bias), f3(r$ate_sd), f3(r$ate_rmse), f2(r$ate_cover), f2(r$ate_power))
    if (sc4) cells <- c(cells, f3(r$surf_rmse_R), f3(r$surf_rmse_Rc), f2(r$map_R), f2(r$map_Rc),
                        f2(r$mapspike_R), f2(r$mapspike_Rc), f2(r$g_mean_R), f2(r$g_mean_Rc),
                        f1(r$ess_prior), f1(r$ess_realized), f1(r$ess_ceiling), f2(r$ess_capped))
    else cells <- c(cells, f2(r$map_all), f2(ifelse(is.na(r$g_mean_R), r$g_mean_Rc, r$g_mean_R)), f1(r$ess_prior), f1(r$ess_realized), f1(r$ess_ceiling), f2(r$ess_capped))
    out <- c(out, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  out <- c(out, "")
}
writeLines(out, file.path(DIR, "tables.md"))

# Go/no-go criteria, spec Section E, LRC-BART-100 unless stated.
g <- function(sc, m, v) full[scenario == sc & method == m][[v]]
crit <- list()
add <- function(id, sc, stmt, value, pass, note = "") crit[[length(crit) + 1]] <<- data.table(id = id, scenario = sc, criterion = stmt, value = value, pass = pass, note = note)
L <- "LRC-BART-100"
for (sc in c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5")) {
  b <- g(sc, L, "ate_bias"); add("Sc3 bias", sc, "|bias| <= 0.05", sprintf("%.3f", b), abs(b) <= 0.05)
  cv <- g(sc, L, "ate_cover"); add("Sc3 coverage", sc, "coverage >= 0.92", sprintf("%.3f", cv), cv >= 0.92)
  r <- g(sc, L, "ate_rmse"); rp <- g(sc, "BART-PP", "ate_rmse")
  add("Sc3 RMSE", sc, "RMSE <= BART-PP + 0.02", sprintf("%.3f vs %.3f + 0.02", r, rp), r <= rp + 0.02)
}
r <- g("Sc1", L, "ate_rmse"); rp <- g("Sc1", "BART-PP", "ate_rmse")
add("Sc1 RMSE", "Sc1", "RMSE <= 0.90 x BART-PP", sprintf("%.3f vs 0.90 x %.3f = %.3f", r, rp, 0.9 * rp), r <= 0.9 * rp)
p <- g("Sc1", L, "ate_power"); pp <- g("Sc1", "BART-PP", "ate_power")
add("Sc1 power", "Sc1", "power >= BART-PP + 0.05", sprintf("%.3f vs %.3f + 0.05", p, pp), p >= pp + 0.05,
    "harness power is rejection of 0 by the 95% interval at delta = 1, near 1 for every method")
sc <- "Sc4_X5_d2"
b <- g(sc, L, "ate_bias"); add("Sc4 d2 bias", sc, "|bias| <= 0.10", sprintf("%.3f", b), abs(b) <= 0.10)
r <- g(sc, L, "surf_rmse_R"); rp <- g(sc, "BART-PP", "surf_rmse_R")
add("Sc4 d2 RMSE_R", sc, "control RMSE in R <= 0.90 x BART-PP's", sprintf("%.3f vs 0.90 x %.3f = %.3f", r, rp, 0.9 * rp), r <= 0.9 * rp)
m <- g(sc, L, "mapspike_R"); add("Sc4 d2 map R (spec)", sc, "all-spike map mean in R >= 0.5", sprintf("%.3f", m), m >= 0.5)
m <- g(sc, L, "mapspike_Rc"); add("Sc4 d2 map Rc (spec)", sc, "all-spike map mean in R^c <= 0.15", sprintf("%.3f", m), m <= 0.15)
m <- g(sc, L, "map_R"); add("Sc4 d2 map R (restated)", sc, "P(|g| < 0.5) mean in R >= 0.7", sprintf("%.3f", m), m >= 0.7)
m <- g(sc, L, "map_Rc"); add("Sc4 d2 map Rc (restated)", sc, "P(|g| < 0.5) mean in R^c <= 0.15", sprintf("%.3f", m), m <= 0.15)
gg <- g(sc, L, "g_mean_R"); add("Sc4 d2 g in R", sc, "mean g within 0.25 of 0 in R", sprintf("%.3f", gg), abs(gg) <= 0.25)
gg <- g(sc, L, "g_mean_Rc"); add("Sc4 d2 g in Rc", sc, "mean g within 0.25 of -2 in R^c", sprintf("%.3f", gg), abs(gg + 2) <= 0.25)
sc <- "Sc4_X5_d1"
b <- g(sc, L, "ate_bias"); add("Sc4 d1 bias", sc, "|bias| <= 0.10", sprintf("%.3f", b), abs(b) <= 0.10)
r <- g(sc, L, "ate_rmse"); rp <- g(sc, "BART-PP", "ate_rmse")
add("Sc4 d1 RMSE", sc, "overall RMSE no worse than BART-PP", sprintf("%.3f vs %.3f", r, rp), r <= rp)
sc <- "Sc5_d2"
b <- g(sc, L, "ate_bias"); add("Sc5 d2 bias", sc, "|bias| <= 0.10", sprintf("%.3f", b), abs(b) <= 0.10)
r <- g(sc, L, "ate_rmse"); rn <- g(sc, "BART-NP", "ate_rmse")
add("Sc5 d2 RMSE", sc, "RMSE <= 1.10 x BART-NP", sprintf("%.3f vs 1.10 x %.3f = %.3f", r, rn, 1.1 * rn), r <= 1.1 * rn)
m <- g(sc, L, "mapspike_Rc"); add("Sc5 d2 map (spec)", sc, "all-spike map mean <= 0.15 (R empty: R^c mean)", sprintf("%.3f", m), m <= 0.15)
m <- g(sc, L, "map_all"); add("Sc5 d2 map (restated)", sc, "P(|g| < 0.5) mean <= 0.15", sprintf("%.3f", m), m <= 0.15)
crit <- rbindlist(crit)
fwrite(crit, file.path(DIR, "criteria.csv"))
cat("| # | scenario | criterion | value | result |\n|---|---|---|---|---|\n")
for (i in seq_len(nrow(crit))) cat(sprintf("| %d | %s | %s | %s | %s |\n", i, crit$scenario[i], crit$criterion[i], crit$value[i], ifelse(crit$pass[i], "PASS", "FAIL")))
cat(sprintf("\nPASS %d of %d\n", sum(crit$pass), nrow(crit)))
