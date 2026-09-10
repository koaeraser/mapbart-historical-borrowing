# make_tables.R: LaTeX tables for Sections 3 and 4 of the LRC-BART revision,
# generated from the Phase 2 full-study summaries (02-validation/full/*/
# summary_table.csv) and the Phase 4 application table
# (04-application/results_table.csv). Every number in the manuscript tables
# comes from these files; nothing is typed by hand. Output: 05-writing/tables/
# tab_*.tex, each a complete table environment to be \input.
#
# Definitions (restated in the text): bias = mean over replicates of the
# posterior mean minus the truth; SD = mean posterior SD; RMSE = Monte Carlo
# RMSE of the posterior mean; Cov = coverage of the central 95% interval;
# Pow = fraction of replicates whose 95% interval excludes the null (0 for the
# ATE, 1 for the survival ratios). Survival point estimates are posterior means.

REV <- Sys.getenv("LRC_REVISION_ROOT", "revision")
OUT <- file.path(REV, "05-writing", "tables")
dir.create(OUT, showWarnings = FALSE)
g  <- read.csv(file.path(REV, "02-validation/full/gaussian/summary_table.csv"), stringsAsFactors = FALSE)
s  <- read.csv(file.path(REV, "02-validation/full/survival/summary_table.csv"), stringsAsFactors = FALSE)
sa <- read.csv(file.path(REV, "02-validation/full/single_arm/summary_table.csv"), stringsAsFactors = FALSE)
ap <- read.csv(file.path(REV, "04-application/results_table.csv"), stringsAsFactors = FALSE)

cahb <- read.csv(file.path(REV, "02-validation/cahb/summary_table.csv"))
regional <- read.csv(file.path(REV, "02-validation/sc4_boundary/regional_ess_summary.csv"))
boundary <- read.csv(file.path(REV, "02-validation/sc4_boundary/summary.csv"))
boundary_paired <- read.csv(file.path(REV, "02-validation/sc4_boundary/paired.csv"))
cahb_rows <- cahb[cahb$method == "CAHB", ]
for (nm in setdiff(names(g), names(cahb_rows))) cahb_rows[[nm]] <- NA
# CAHB diagnostics have different definitions; do not print them as B_c or LRC ESS.
g_main <- rbind(g, cahb_rows[, names(g)])
g_main$map_R[g_main$method == "CAHB"] <- NA
g_main$map_Rc[g_main$method == "CAHB"] <- NA
g_main$ess_realized[g_main$method == "CAHB"] <- NA

# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "---", sub("^-", "$-$", formatC(x, format = "f", digits = 3)))
f2 <- function(x) ifelse(is.na(x), "---", sub("^-", "$-$", formatC(x, format = "f", digits = 2)))
f1 <- function(x) ifelse(is.na(x), "---", formatC(x, format = "f", digits = 1))
f0 <- function(x) ifelse(is.na(x), "---", formatC(x, format = "f", digits = 0))
# very large values (the hierarchical single-arm rows) are printed as ">10^3"
fbig <- function(x, fmt = f3) ifelse(!is.na(x) & abs(x) > 1000, "$>\\!10^{3}$", fmt(x))

method_name <- c(
  "CAHB" = "CAHB", "CAHB-p10" = "CAHB-p10", "CAHB-lit" = "CAHB-lit",
  "LM-NP" = "LM-NP", "LM-CP" = "LM-CP", "LM-PP" = "LM-PP",
  "AFT-NP" = "AFT-NP", "AFT-CP" = "AFT-CP", "AFT-PP" = "AFT-PP",
  "BART-NP" = "BART-NP", "BART-CP" = "BART-CP", "BART-PP" = "BART-PP",
  "PSCL" = "PSCL", "MAP-100" = "MAP (100)", "MAP-75" = "MAP (75)", "MAP-50" = "MAP (50)",
  "HierLM" = "HierLM", "HierAFT" = "HierAFT",
  "C-BART" = "C-BART", "LRC-BART-w0" = "LRC-BART ($w=0$)",
  "LRC-BART-100" = "LRC-BART (100)", "LRC-BART-75" = "LRC-BART (75)", "LRC-BART-50" = "LRC-BART (50)",
  "LRC-BART-f90" = "LRC-BART (0.9)", "LRC-BART-f50" = "LRC-BART (0.5)", "LRC-BART-f25" = "LRC-BART (0.25)",
  "LRC-BART-Hg10" = "LRC-BART ($H_g=10$)",
  "LRC-BART-full" = "LRC-BART (pooled)", "LRC-BART-100-w1" = "LRC-BART (100)",
  "LRC-BART-f90-w1" = "LRC-BART (0.9)", "LRC-BART-100-w0.9" = "LRC-BART ($w=0.9$)",
  "LRC-BART-f90-w0.9" = "LRC-BART (0.9, $w=0.9$)")
scen_name <- c(
  "Sc1" = "Sc1", "Sc2" = "Sc2",
  "Sc3_rho-0.5" = "Sc3, $\\rho=-0.5$", "Sc3_rho0" = "Sc3, $\\rho=0$", "Sc3_rho0.5" = "Sc3, $\\rho=0.5$",
  "Sc4_X5_d0.5" = "Sc4, $X_5$, $\\delta_2=0.5$", "Sc4_X5_d1" = "Sc4, $X_5$, $\\delta_2=1$",
  "Sc4_X5_d2" = "Sc4, $X_5$, $\\delta_2=2$", "Sc4_X7_d1" = "Sc4, $X_7$, $\\delta_2=1$",
  "Sc4_X7_d2" = "Sc4, $X_7$, $\\delta_2=2$", "Sc5_d1" = "Sc5, $\\delta_2=1$", "Sc5_d2" = "Sc5, $\\delta_2=2$")
scen_short <- c(
  "Sc1" = "Sc1", "Sc2" = "Sc2", "Sc3_rho-0.5" = "Sc3 ($-0.5$)", "Sc3_rho0" = "Sc3 (0)", "Sc3_rho0.5" = "Sc3 (0.5)",
  "Sc4_X5_d0.5" = "Sc4 (0.5)", "Sc4_X5_d1" = "Sc4 (1)", "Sc4_X5_d2" = "Sc4 (2)",
  "Sc4_X7_d1" = "Sc4$'$ (1)", "Sc4_X7_d2" = "Sc4$'$ (2)", "Sc5_d1" = "Sc5 (1)", "Sc5_d2" = "Sc5 (2)")

wrap_table <- function(body, caption, label, colspec, header, size = "\\scriptsize", stretch = 0.8,
                       colsep = "4pt", env = "table", note = NULL, resize = FALSE) {
  if (label %in% c("tab:gauss", "tab:map", "tab:single")) {
    # Keep complete scenario groups together and let the table span pages.
    isrow <- grepl("\\\\\\\\$", body)
    for (i in which(isrow)) if (i < length(body) && body[i + 1] != "\\midrule") body[i] <- paste0(body[i], "*")
    return(c("\\begingroup", "\\singlespacing", "\\scriptsize",
      "\\setlength{\\tabcolsep}{2.5pt}", "\\renewcommand{\\arraystretch}{1.0}",
      sprintf("\\begin{longtable}{%s}", colspec),
      sprintf("\\caption{%s}\\label{%s}\\\\", caption, label),
      "\\toprule", header, "\\midrule", "\\endfirsthead",
      sprintf("\\multicolumn{%d}{l}{\\textit{Table \\ref{%s}, continued}}\\\\", nchar(colspec), label),
      "\\toprule", header, "\\midrule", "\\endhead",
      "\\midrule", "\\endfoot", "\\bottomrule", "\\endlastfoot",
      body, "\\end{longtable}", "\\endgroup"))
  }
  c(sprintf("\\begin{%s}[htbp]", env), "\\centering", "\\singlespacing",
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label), size,
    sprintf("\\setlength{\\tabcolsep}{%s}", colsep),
    sprintf("\\renewcommand{\\arraystretch}{%s}", stretch),
    sprintf("\\begin{adjustbox}{max width=\\textwidth,max totalheight=%s\\textheight}", if (label == "tab:map") "0.50" else "0.68"),
    sprintf("\\begin{tabular}{%s}", colspec), "\\toprule", header, "\\midrule",
    body, "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}",
    if (!is.null(note)) c("\\par\\vspace{2pt}", sprintf("\\parbox{\\textwidth}{\\scriptsize %s}", note)),
    sprintf("\\end{%s}", env))
}
block_rows <- function(df, scen_col, row_fun, scen_labels, scen_order, method_order) {
  out <- character()
  first <- TRUE
  for (sc in scen_order) {
    d <- df[df[[scen_col]] == sc, ]
    d <- d[match(method_order, d$method), ]
    d <- d[!is.na(d$method), ]
    if (nrow(d) == 0) next
    if (!first) out <- c(out, "\\midrule")
    first <- FALSE
    rows <- row_fun(d)
    rows[1] <- sprintf("\\multirow{%d}{*}{%s} %s", nrow(d), scen_labels[sc], rows[1])
    out <- c(out, rows)
  }
  out
}

# ---------------------------------------------------------------------------
# Table: Gaussian two-arm, main text (selected scenarios and methods)
# ---------------------------------------------------------------------------
g_main_sc <- c("Sc1", "Sc2", "Sc3_rho0", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d2")
g_main_m  <- c("LM-PP", "BART-NP", "BART-CP", "BART-PP", "PSCL", "MAP-100", "HierLM", "C-BART",
               "LRC-BART-100", "LRC-BART-f90", "CAHB")
row_g <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], f3(d$bias), f3(d$sd), f3(d$rmse), f2(d$cov), f2(d$pow),
  f3(d$ctrl_bias), f3(d$ctrl_rmse))
hdr_g <- c("& & \\multicolumn{5}{c}{ATE} & \\multicolumn{2}{c}{Control mean} \\\\",
           "\\cmidrule(lr){3-7}\\cmidrule(lr){8-9}",
           "Scenario & Method & Bias & SD & RMSE & Cov & Pow & Bias & RMSE \\\\")
cap_g <- paste0("Gaussian outcomes, two-arm design, $R=500$ replicates: average treatment effect (ATE) ",
  "bias, mean posterior SD, Monte Carlo RMSE, coverage of the 95\\% credible interval, and power ",
  "(95\\% interval excludes zero at $\\delta=1$), with bias and RMSE of the RCT-standardized control mean. ",
  "Sc4 raises the RWD control surface by $\\delta_2$ outside $\\mathcal R=\\{X_5>2\\}$; Sc5 raises it everywhere. ",
  "LRC-BART (100) targets a prior ESS of 100 RCT controls and LRC-BART (0.9) targets 0.9 of the ESS ceiling; ",
  "CAHB uses $R=200$ paired replicates and our transfer described in Web Appendix~G.2; other rows use $R=500$. C-BART is LRC-BART with $w=1$. PSCL reports no arm-specific means. The full table with every scenario and method is Web Table~S1.")
writeLines(wrap_table(block_rows(g_main, "scenario", row_g, scen_name, g_main_sc, g_main_m), cap_g,
                      "tab:gauss", "llccccccc", hdr_g, colsep = "3pt"), file.path(OUT, "tab_gauss.tex"))

# full Gaussian table (Web Appendix), all scenarios, all methods, landscape-free with longtable
g_all_m <- c("LM-NP", "LM-CP", "LM-PP", "BART-NP", "BART-CP", "BART-PP", "PSCL", "MAP-100", "MAP-75", "MAP-50",
             "HierLM", "C-BART", "LRC-BART-w0", "LRC-BART-100", "LRC-BART-75", "LRC-BART-50",
             "LRC-BART-f90", "LRC-BART-f50", "LRC-BART-f25", "LRC-BART-Hg10")
g_all_sc <- names(scen_name)
row_gf <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], f3(d$bias), f3(d$sd), f3(d$rmse), f2(d$cov), f2(d$pow),
  f3(d$trt_bias), f3(d$trt_rmse), f3(d$ctrl_bias), f3(d$ctrl_rmse))
hdr_gf <- c("& & \\multicolumn{5}{c}{ATE} & \\multicolumn{2}{c}{Treated mean} & \\multicolumn{2}{c}{Control mean} \\\\",
            "\\cmidrule(lr){3-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}",
            "Scenario & Method & Bias & SD & RMSE & Cov & Pow & Bias & RMSE & Bias & RMSE \\\\")
body_gf <- block_rows(g, "scenario", row_gf, scen_name, g_all_sc, g_all_m)
cap_gf <- paste0("Gaussian outcomes, two-arm design, all scenarios and methods, $R=500$. Columns as in Table~\\ref{tab:gauss}, ",
  "plus bias and RMSE of the RCT-standardized treated mean. MAP ($N$) is the meta-analytic predictive prior calibrated to prior ESS $N$; ",
  "LRC-BART ($N$) targets prior ESS $N$ and LRC-BART ($p$) the fraction $p$ of the ESS ceiling; LRC-BART ($w=0$) fixes every leaf in the slab; ",
  "LRC-BART ($H_g=10$) doubles the discrepancy ensemble.")
writeLines(c("\\begin{center}", "\\scriptsize", "\\setlength{\\tabcolsep}{3.5pt}", "\\renewcommand{\\arraystretch}{0.8}",
             "\\begin{longtable}{llccccccccc}",
             sprintf("\\caption{%s}\\label{tab:gauss_full}\\\\", cap_gf),
             "\\toprule", hdr_gf, "\\midrule", "\\endfirsthead",
             "\\multicolumn{11}{l}{\\textit{Web Table S1, continued}}\\\\", "\\toprule", hdr_gf, "\\midrule", "\\endhead",
             "\\midrule", "\\multicolumn{11}{r}{\\textit{continued on the next page}}\\\\", "\\endfoot",
             "\\bottomrule", "\\endlastfoot",
             body_gf, "\\end{longtable}", "\\end{center}"), file.path(OUT, "tab_gauss_full.tex"))

# ---------------------------------------------------------------------------
# Table: Sc4 and Sc5 borrowing map and control surface (Gaussian)
# ---------------------------------------------------------------------------
m_sc <- c("Sc1", "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc4_X7_d1", "Sc4_X7_d2", "Sc5_d1", "Sc5_d2")
m_m  <- c("BART-PP", "C-BART", "LRC-BART-w0", "LRC-BART-100", "CAHB")
row_m <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], f3(d$bias), f2(d$cov),
  f3(d$surf_rmse_R), f3(d$surf_rmse_Rc),
  f2(d$map_R), f2(d$map_Rc), f2(d$g_mean_R), f2(d$g_mean_Rc), f1(d$ess_realized))
hdr_m <- c("& & \\multicolumn{2}{c}{ATE} & \\multicolumn{2}{c}{Surface RMSE} & \\multicolumn{2}{c}{$\\bar{\\mathcal B}_c$} & \\multicolumn{2}{c}{Mean $\\hat g$} & \\\\",
           "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}",
           "Scenario & Method & Bias & Cov & $\\mathcal R$ & $\\mathcal R^c$ & $\\mathcal R$ & $\\mathcal R^c$ & $\\mathcal R$ & $\\mathcal R^c$ & ESS \\\\")
cap_m <- paste0("Regional and global shifts, Gaussian outcomes, $R=500$: ATE bias and coverage; RMSE of the posterior mean control surface ",
  "over all $N$ RCT profiles (treated and control) inside the compatible region $\\mathcal R$ and outside it; the borrowing map ",
  "$\\mathcal B_c(x)=P(|g(x)|<0.5\\mid\\text{data})$ averaged over the RCT control profiles in each region; the posterior mean of $g(x)$ ",
  "averaged over each region (true value $0$ in $\\mathcal R$ and $-\\delta_2$ in $\\mathcal R^c$); and the realized ESS. ",
  "In Sc5 the compatible region is empty. Sc4$'$ places the shift on the null covariate $X_7$. LRC-BART ($w=0$) has no spike, so its map ",
  "reflects the slab alone. CAHB rows use $R=200$; all other rows use $R=500$. CAHB precision shares and effective sizes have different definitions and are reported separately in Web Appendix~G.2. For Sc1, the two columns use the Sc4 partition although both regions are compatible.")
writeLines(wrap_table(block_rows(g_main, "scenario", row_m, scen_short, m_sc, m_m), cap_m,
                      "tab:map", "llccccccccc", hdr_m, colsep = "3pt", resize = TRUE), file.path(OUT, "tab_map.tex"))
mapfile <- file.path(OUT, "tab_map.tex")
maplines <- readLines(mapfile)
reg_sc <- c("Sc1", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d2")
regrows <- vapply(reg_sc, function(sc) {
  d <- regional[regional$scenario == sc, ]
  sprintf("%s & %s & %s \\\\", scen_short[sc],
          f2(d$ess_real_S[d$region == "R"]), f2(d$ess_real_S[d$region == "Rc"]))
}, "")
regpanel <- c("\\par\\medskip", "\\begin{tabular}{lcc}", "\\toprule",
  "Scenario & Regional ESS, $X_5>2$ & Regional ESS, $X_5\\le2$ \\\\", "\\midrule",
  regrows, "\\bottomrule", "\\end{tabular}",
  "\\par\\smallskip\\parbox{\\textwidth}{\\scriptsize Regional realized ESS for LRC-BART (100), averaged over the first 200 paired replicates. Each estimand is the control mean standardized to all RCT profiles in that region. The partition is held fixed for Sc1 and Sc5; both regions are compatible in Sc1 and neither is compatible in Sc5. Regional effective sizes are not additive.}")
writeLines(c(maplines, "\\begin{center}\\scriptsize", regpanel, "\\end{center}"), mapfile)


# ---------------------------------------------------------------------------
# Table: fractional ladder (Gaussian, Web Appendix)
# ---------------------------------------------------------------------------
l_m <- c("LRC-BART-50", "LRC-BART-75", "LRC-BART-100", "LRC-BART-f25", "LRC-BART-f50", "LRC-BART-f90")
row_l <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], f3(d$bias), f3(d$rmse), f2(d$cov), f1(d$ess_prior), f1(d$ess_realized),
  f0(d$ess_ceiling), f2(d$ess_capped))
hdr_l <- "Scenario & Target & Bias & RMSE & Cov & Prior ESS & Realized ESS & Ceiling & Capped \\\\"
cap_l <- paste0("The ESS ladder, Gaussian outcomes, $R=500$: absolute targets $N_{\\rm target}\\in\\{50,75,100\\}$ and fractional ",
  "targets $0.25$, $0.5$, $0.9$ of the ceiling $\\sigma_1^2/V_\\mu^f$. Prior ESS is $\\mathrm{ESS}_0$ at the calibrated scale, ",
  "realized ESS its posterior counterpart, ceiling the mean ceiling, and capped the fraction of replicates in which the target ",
  "was at or above the ceiling and the calibration returned the point-mass limit. Prior ESS is a block average and can exceed the whole-chain ceiling by a few percent when capped (Section~\\ref{sec:sim_design}).")
writeLines(wrap_table(block_rows(g, "scenario", row_l, scen_name, c("Sc1", "Sc2", "Sc3_rho0", "Sc4_X5_d1", "Sc4_X5_d2", "Sc5_d2"), l_m),
                      cap_l, "tab:ladder", "llccccccc", hdr_l), file.path(OUT, "tab_ladder.tex"))

# ---------------------------------------------------------------------------
# Table: survival two-arm, main text (median ratio and RMST ratio side by side)
# ---------------------------------------------------------------------------
s_med <- s[s$estimand == "med", ]; s_rm <- s[s$estimand == "rmst", ]
s_main_sc <- c("Sc1", "Sc2", "Sc3_rho0", "Sc4_X5_d1", "Sc5_d1")
s_main_m  <- c("AFT-PP", "BART-NP", "BART-CP", "BART-PP", "HierAFT", "C-BART", "LRC-BART-100", "LRC-BART-f90")
row_s <- function(d) {
  r <- s_rm[match(paste(d$scenario, d$method), paste(s_rm$scenario, s_rm$method)), ]
  sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          method_name[d$method], f3(d$bias), f3(d$rmse), f2(d$cov), f2(d$pow),
          f3(r$bias), f3(r$rmse), f2(r$cov), f2(r$pow))
}
hdr_s <- c("& & \\multicolumn{4}{c}{Median ratio} & \\multicolumn{4}{c}{RMST ratio, $t^\\ast=3$} \\\\",
           "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
           "Scenario & Method & Bias & RMSE & Cov & Pow & Bias & RMSE & Cov & Pow \\\\")
cap_s <- paste0("Survival outcomes, two-arm design, $R=500$: bias, Monte Carlo RMSE, coverage and power (95\\% interval excludes one) ",
  "for the ratio of standardized population medians (true value $1.4$) and the ratio of standardized restricted mean survival times at ",
  "$t^\\ast=3$. Sc4 raises the RWD log-time surface by $\\delta_2=1$ outside $\\mathcal R=\\{X_5>2\\}$ and Sc5 by $1$ everywhere. ",
  "Point estimates are posterior means. The full table is Web Table~S2.")
writeLines(wrap_table(block_rows(s_med, "scenario", row_s, scen_name, s_main_sc, s_main_m), cap_s,
                      "tab:surv", "llcccccccc", hdr_s), file.path(OUT, "tab_surv.tex"))

s_all_m <- c("AFT-NP", "AFT-CP", "AFT-PP", "BART-NP", "BART-CP", "BART-PP", "HierAFT", "C-BART", "LRC-BART-w0",
             "LRC-BART-100", "LRC-BART-75", "LRC-BART-50", "LRC-BART-f90", "LRC-BART-f50", "LRC-BART-f25")
s_all_sc <- c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5", "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc5_d1")
row_sf <- function(d) {
  r <- s_rm[match(paste(d$scenario, d$method), paste(s_rm$scenario, s_rm$method)), ]
  sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          method_name[d$method], f3(d$bias), f3(d$sd), f3(d$rmse), f2(d$cov), f2(d$pow), f3(d$ctrl_bias),
          f3(r$bias), f3(r$sd), f3(r$rmse), f2(r$cov), f2(r$pow))
}
hdr_sf <- c("& & \\multicolumn{6}{c}{Median ratio} & \\multicolumn{5}{c}{RMST ratio, $t^\\ast=3$} \\\\",
            "\\cmidrule(lr){3-8}\\cmidrule(lr){9-13}",
            "Scenario & Method & Bias & SD & RMSE & Cov & Pow & Ctrl bias & Bias & SD & RMSE & Cov & Pow \\\\")
cap_sf <- paste0("Survival outcomes, two-arm design, all scenarios and methods, $R=500$. Columns as in Table~\\ref{tab:surv}, with the mean ",
  "posterior SD and the bias of the standardized control median on the log scale (Ctrl bias).")
writeLines(c("\\begin{center}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\renewcommand{\\arraystretch}{0.8}",
             "\\begin{longtable}{llccccccccccc}",
             sprintf("\\caption{%s}\\label{tab:surv_full}\\\\", cap_sf),
             "\\toprule", hdr_sf, "\\midrule", "\\endfirsthead",
             "\\multicolumn{13}{l}{\\textit{Web Table S2, continued}}\\\\", "\\toprule", hdr_sf, "\\midrule", "\\endhead",
             "\\midrule", "\\multicolumn{13}{r}{\\textit{continued on the next page}}\\\\", "\\endfoot",
             "\\bottomrule", "\\endlastfoot",
             block_rows(s_med, "scenario", row_sf, scen_name, s_all_sc, s_all_m),
             "\\end{longtable}", "\\end{center}"), file.path(OUT, "tab_surv_full.tex"))

# survival map table (Web Appendix)
sm_sc <- c("Sc4_X5_d0.5", "Sc4_X5_d1", "Sc5_d1"); sm_m <- c("BART-PP", "C-BART", "LRC-BART-w0", "LRC-BART-100")
row_sm <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], f3(d$bias), f2(d$cov), f3(d$surf_rmse_R), f3(d$surf_rmse_Rc),
  f2(d$map_R), f2(d$map_Rc), f2(d$g_mean_R), f2(d$g_mean_Rc), f1(d$ess_realized))
cap_sm <- paste0("Regional and global shifts, survival outcomes (median ratio), $R=500$: columns as in Table~\\ref{tab:map}, ",
  "with the surface RMSE, the map threshold $c=0.5$ and $g$ on the log-time scale.")
writeLines(wrap_table(block_rows(s_med, "scenario", row_sm, scen_short, sm_sc, sm_m), cap_sm,
                      "tab:map_surv", "llccccccccc", hdr_m, colsep = "3pt", resize = TRUE), file.path(OUT, "tab_map_surv.tex"))

# ---------------------------------------------------------------------------
# Table: single-arm design
# ---------------------------------------------------------------------------
sa$key <- paste(sa$design, sa$estimand)
sa_blocks <- c("gaussian_n200 ate" = "Gaussian, $n_T=200$, ATE",
               "survival_n30 med" = "Survival, $n_T=30$, median ratio",
               "survival_n30 rmst" = "Survival, $n_T=30$, RMST ratio",
               "survival_n200 med" = "Survival, $n_T=200$, median ratio",
               "survival_n200 rmst" = "Survival, $n_T=200$, RMST ratio")
sa_full_m <- c("LM-CP", "AFT-CP", "HierLM", "HierAFT", "BART-CP", "LRC-BART-full", "LRC-BART-100-w1", "LRC-BART-f90-w1", "LRC-BART-100-w0.9", "LRC-BART-f90-w0.9")
sa_main_m <- c("LM-CP", "AFT-CP", "HierLM", "HierAFT", "BART-CP", "LRC-BART-full", "LRC-BART-f90-w1", "LRC-BART-100-w0.9")
sa_blocks2 <- c("gaussian_n200 ate" = "Gaussian, $n_T=200$\\\\ATE",
                "survival_n30 med" = "Survival, $n_T=30$\\\\median ratio",
                "survival_n30 rmst" = "Survival, $n_T=30$\\\\RMST ratio",
                "survival_n200 med" = "Survival, $n_T=200$\\\\median ratio",
                "survival_n200 rmst" = "Survival, $n_T=200$\\\\RMST ratio")
row_sa <- function(d) sprintf("& %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  method_name[d$method], fbig(d$bias), fbig(d$sd), fbig(d$rmse), f2(d$cov), f2(d$pow),
  fbig(d$ctrl_bias), fbig(d$ctrl_rmse), f1(d$ess_prior))
sa_body <- function(methods, blocks = names(sa_blocks2)) {
  out <- character(); first <- TRUE
  for (k in blocks) {
    for (sc in c("Sc1", "Sc2")) {
      d <- sa[sa$key == k & sa$scenario == sc, ]
      d <- d[match(methods, d$method), ]; d <- d[!is.na(d$method), ]
      if (!first) out <- c(out, "\\midrule"); first <- FALSE
      rows <- row_sa(d)
      rows[1] <- sprintf("\\multirow{%d}{*}{\\shortstack[l]{%s\\\\%s}} %s", nrow(d), sa_blocks2[k], sc, rows[1])
      out <- c(out, rows)
    }
  }
  out
}
hdr_sa <- c("& & \\multicolumn{5}{c}{Treatment effect} & \\multicolumn{2}{c}{Control} & \\\\",
            "\\cmidrule(lr){3-7}\\cmidrule(lr){8-9}",
            "Design & Method & Bias & SD & RMSE & Cov & Pow & Bias & RMSE & Prior ESS \\\\")
cap_sa <- paste0("Single-arm design, $R=500$: the trial contributes $n_T$ treated subjects only and the $n_2=300$ RWD controls are the sole ",
  "control source, under Sc1 (no confounding) and Sc2 (measured confounding). LRC-BART (pooled) is $w=1$ with the spike scale at the ",
  "grid minimum, LRC-BART (0.9) is $w=1$ at 0.9 of the ESS ceiling, and LRC-BART ($w=0.9$) is the slab sensitivity at the absolute target 100. ",
  "HierLM and HierAFT use the stated hyperprior $\\tau_j^2\\sim$ Inv-Gamma$(1.5,0.75)$ with no RCT controls, under which the control ",
  "surface is drawn from the hierarchy; values above $10^3$ are printed as such. Control columns are on the outcome scale (Gaussian) or ",
  "the log-median scale (survival). Web Table~S6 adds the absolute-target rows and the median ratio at $n_T=200$.")
writeLines(wrap_table(sa_body(sa_main_m, setdiff(names(sa_blocks2), "survival_n200 med")), cap_sa, "tab:single", "llcccccccc", hdr_sa, stretch = 0.8, colsep = "3pt", resize = TRUE),
           file.path(OUT, "tab_single.tex"))
cap_saf <- paste0("Single-arm design, all methods, $R=500$. Columns as in Table~\\ref{tab:single}; LRC-BART (100) and LRC-BART ($w=0.9$) ",
  "use the absolute target 100, capped at the ceiling where infeasible, and LRC-BART (0.9, $w=0.9$) the fractional target with the slab.")
writeLines(c("\\begin{center}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\renewcommand{\\arraystretch}{0.8}",
             "\\begin{longtable}{llcccccccc}",
             sprintf("\\caption{%s}\\label{tab:single_full}\\\\", cap_saf),
             "\\toprule", hdr_sa, "\\midrule", "\\endfirsthead",
             "\\multicolumn{10}{l}{\\textit{Web Table S6, continued}}\\\\", "\\toprule", hdr_sa, "\\midrule", "\\endhead",
             "\\midrule", "\\multicolumn{10}{r}{\\textit{continued on the next page}}\\\\", "\\endfoot",
             "\\bottomrule", "\\endlastfoot",
             sa_body(sa_full_m), "\\end{longtable}", "\\end{center}"), file.path(OUT, "tab_single_full.tex"))

# ---------------------------------------------------------------------------
# Table: application
# ---------------------------------------------------------------------------
ci <- function(m, lo, hi, med = NULL) {
  if (is.null(med)) sprintf("%s (%s, %s)", f2(m), f2(lo), f2(hi))
  else sprintf("%s [%s] (%s, %s)", f2(m), f2(med), f2(lo), f2(hi))
}
pick <- function(method, source, coding, est) {
  r <- ap[ap$method == method & grepl(source, ap$source, fixed = TRUE) & ap$coding == coding & ap$estimand == est, ]
  stopifnot(nrow(r) == 1); r
}
rows <- character()
add <- function(label, r, with_med = TRUE, ess = "") {
  rows <<- c(rows, sprintf("%s & %s & %s & %s \\\\", label,
    if (is.na(r$pfs_mean)) ci(r$pfs_median, r$pfs_lo, r$pfs_hi) else ci(r$pfs_mean, r$pfs_lo, r$pfs_hi, r$pfs_median),
    if (is.na(r$os_mean)) ci(r$os_median, r$os_lo, r$os_hi) else ci(r$os_mean, r$os_lo, r$os_hi, r$os_median), ess))
}
rr <- "rerun"
rows <- c(rows, "\\multicolumn{4}{l}{\\textit{Reference methods, harmonized coding}} \\\\")
add("Kaplan--Meier, unadjusted", pick("KM (unadjusted)", rr, "none", "rmst"))
add("AFT, complete pooling", pick("AFT, complete pooling", rr, "harm", "rmst"))
add("AFT, hierarchical ($s^2=0.5$)", pick("AFT, hierarchical (s2 = 0.5)", rr, "harm", "rmst"))
add("AFT-BART on UCMM, no discrepancy", pick("Standard BART (lrcbart plain AFT-BART, 50 trees)", rr, "harm", "rmst"))
rows <- c(rows, "\\midrule", "\\multicolumn{4}{l}{\\textit{LRC-BART, single-arm mode, harmonized coding, $H_f=10$, $w=1$}} \\\\")
for (tg in c("N_target = 30", "N_target = 22.5", "N_target = 15")) {
  r <- pick(sprintf("LRC-BART single-arm (%s, H_f = 10, H_g = 5, w = 1)", tg), "this rerun", "harm", "rmst")
  add(sprintf("$N_{\\rm target}=%s$", sub("N_target = ", "", tg)), r,
      ess = sprintf("%s / %s", f1(r$ess0_pfs), f1(r$ess0_os)))
}
for (tg in c("0.9 x ceiling", "0.5 x ceiling", "0.25 x ceiling")) {
  r <- pick(sprintf("LRC-BART single-arm (%s, H_f = 10, H_g = 5, w = 1)", tg), "this rerun", "harm", "rmst")
  add(sprintf("$%s\\times$ ceiling", sub(" x ceiling", "", tg)), r,
      ess = sprintf("%s / %s", f1(r$ess0_pfs), f1(r$ess0_os)))
}
rows <- c(rows, "\\midrule", "\\multicolumn{4}{l}{\\textit{Sensitivity analyses at $N_{\\rm target}=30$}} \\\\")
r <- pick("LRC-BART single-arm (N_target = 30, H_f = 10, H_g = 5, w = 0.9)", "this rerun", "harm", "rmst")
add("$w=0.9$", r, ess = sprintf("%s / %s", f1(r$ess0_pfs), f1(r$ess0_os)))
r <- pick("LRC-BART single-arm (N_target = 30, H_f = 50, H_g = 5, w = 1)", "this rerun", "harm", "rmst")
add("$H_f=50$", r, ess = sprintf("%s / %s", f1(r$ess0_pfs), f1(r$ess0_os)))
r <- pick("LRC-BART single-arm (N_target = 30, H_f = 10, H_g = 5, w = 1)", "this rerun", "orig", "rmst")
add("Original covariate coding", r, ess = sprintf("%s / %s", f1(r$ess0_pfs), f1(r$ess0_os)))
rows <- c(rows, "\\midrule", "\\multicolumn{4}{l}{\\textit{Published rows (original coding, posterior medians)}} \\\\")
pp <- "paper"
add("AFT, complete pooling", pick("AFT, complete pooling", pp, "orig", "rmst"))
add("AFT, hierarchical ($s^2=0.5$)", pick("AFT, hierarchical (s2 = 0.5)", pp, "orig", "rmst"))
add("Standard BART", pick("Standard BART", pp, "orig", "rmst"))
add("MAP-AFT-BART, $N_{\\rm target}=30$", pick("MAP-AFT-BART (N_target = 30)", pp, "orig", "rmst"))
add("MAP-AFT-BART, $N_{\\rm target}=15$", pick("MAP-AFT-BART (N_target = 15)", pp, "orig", "rmst"))
hdr_a <- "Method & PFS & OS & $\\mathrm{ESS}_0$ (PFS / OS) \\\\"
cap_a <- paste0("EloKRd versus UCMM: standardized 5-year RMST ratio, posterior mean [posterior median] and 95\\% credible interval, ",
  "for progression-free survival (PFS) and overall survival (OS) at the 30 EloKRd profiles. The ESS ceilings are 63.7 (PFS) and 36.7 (OS). ",
  "Reference rows are rerun under the harmonized covariate coding of Section~\\ref{sec:app_data}; the Kaplan--Meier row is unadjusted ",
  "and reports the estimate with its large-sample interval from the survRM2 package. The published rows are those of the earlier analysis under the original coding, ",
  "reported there as posterior medians.")
writeLines(wrap_table(rows, cap_a, "tab:app", "lccc", hdr_a, size = "\\small", stretch = 1.0, colsep = "5pt", resize = TRUE),
           file.path(OUT, "tab_app.tex"))

# ---------------------------------------------------------------------------
# numbers quoted in the text, for the consistency audit
# ---------------------------------------------------------------------------
q <- function(df, sc, m, col, est = NULL) {
  d <- df[df$scenario == sc & df$method == m, ]; if (!is.null(est)) d <- d[d$estimand == est, ]
  d[[col]]
}
quotes <- list(
  gauss_Sc1_rmse_LRC100 = q(g, "Sc1", "LRC-BART-100", "rmse"), gauss_Sc1_rmse_BARTPP = q(g, "Sc1", "BART-PP", "rmse"),
  gauss_Sc1_rmse_BARTNP = q(g, "Sc1", "BART-NP", "rmse"), gauss_Sc1_rmse_BARTCP = q(g, "Sc1", "BART-CP", "rmse"),
  gauss_Sc1_rmse_CBART = q(g, "Sc1", "C-BART", "rmse"),
  gauss_Sc1_ess_prior_LRC100 = q(g, "Sc1", "LRC-BART-100", "ess_prior"), gauss_Sc1_ess_real_LRC100 = q(g, "Sc1", "LRC-BART-100", "ess_realized"),
  gauss_Sc1_ceiling = q(g, "Sc1", "LRC-BART-100", "ess_ceiling"), gauss_Sc1_capped = q(g, "Sc1", "LRC-BART-100", "ess_capped"),
  gauss_Sc2_ceiling = q(g, "Sc2", "LRC-BART-100", "ess_ceiling"),
  gauss_Sc2_bias_LRC100 = q(g, "Sc2", "LRC-BART-100", "bias"), gauss_Sc2_bias_BARTCP = q(g, "Sc2", "BART-CP", "bias"),
  gauss_Sc2_bias_CBART = q(g, "Sc2", "C-BART", "bias"), gauss_Sc2_bias_f25 = q(g, "Sc2", "LRC-BART-f25", "bias"),
  gauss_Sc3_bias_LRC100 = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(g, sc, "LRC-BART-100", "bias")),
  gauss_Sc3_bias_BARTPP = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(g, sc, "BART-PP", "bias")),
  gauss_Sc3_rmse_LRC100 = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(g, sc, "LRC-BART-100", "rmse")),
  gauss_Sc3_rmse_BARTPP = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(g, sc, "BART-PP", "rmse")),
  gauss_Sc3_ess_real_LRC100 = q(g, "Sc3_rho0", "LRC-BART-100", "ess_realized"),
  gauss_Sc3_map_LRC100 = q(g, "Sc3_rho0", "LRC-BART-100", "map_all"),
  gauss_Sc3_bias_CBART = q(g, "Sc3_rho0", "C-BART", "bias"), gauss_Sc3_bias_BARTCP = q(g, "Sc3_rho0", "BART-CP", "bias"),
  gauss_Sc3_bias_Hg10 = q(g, "Sc3_rho0", "LRC-BART-Hg10", "bias"),
  gauss_Sc4d2_map = c(q(g, "Sc4_X5_d2", "LRC-BART-100", "map_R"), q(g, "Sc4_X5_d2", "LRC-BART-100", "map_Rc")),
  gauss_Sc4d1_map = c(q(g, "Sc4_X5_d1", "LRC-BART-100", "map_R"), q(g, "Sc4_X5_d1", "LRC-BART-100", "map_Rc")),
  gauss_Sc4d05_map = c(q(g, "Sc4_X5_d0.5", "LRC-BART-100", "map_R"), q(g, "Sc4_X5_d0.5", "LRC-BART-100", "map_Rc")),
  gauss_Sc4d2_g = c(q(g, "Sc4_X5_d2", "LRC-BART-100", "g_mean_R"), q(g, "Sc4_X5_d2", "LRC-BART-100", "g_mean_Rc")),
  gauss_Sc4d2_surf_LRC = c(q(g, "Sc4_X5_d2", "LRC-BART-100", "surf_rmse_R"), q(g, "Sc4_X5_d2", "LRC-BART-100", "surf_rmse_Rc")),
  gauss_Sc4d2_surf_PP = c(q(g, "Sc4_X5_d2", "BART-PP", "surf_rmse_R"), q(g, "Sc4_X5_d2", "BART-PP", "surf_rmse_Rc")),
  gauss_Sc4d2_surf_CBART_Rc = q(g, "Sc4_X5_d2", "C-BART", "surf_rmse_Rc"),
  gauss_Sc4d1_bias_LRC = q(g, "Sc4_X5_d1", "LRC-BART-100", "bias"), gauss_Sc4d1_rmse_LRC = q(g, "Sc4_X5_d1", "LRC-BART-100", "rmse"),
  gauss_Sc4d1_rmse_PP = q(g, "Sc4_X5_d1", "BART-PP", "rmse"), gauss_Sc4d1_cov_LRC = q(g, "Sc4_X5_d1", "LRC-BART-100", "cov"),
  gauss_Sc4d1_bias_CBART = q(g, "Sc4_X5_d1", "C-BART", "bias"), gauss_Sc4d1_bias_BARTCP = q(g, "Sc4_X5_d1", "BART-CP", "bias"),
  gauss_Sc4d1_bias_f25 = q(g, "Sc4_X5_d1", "LRC-BART-f25", "bias"),
  gauss_Sc4d2_bias_LRC = q(g, "Sc4_X5_d2", "LRC-BART-100", "bias"), gauss_Sc4d2_bias_CBART = q(g, "Sc4_X5_d2", "C-BART", "bias"),
  gauss_Sc4d2_bias_BARTCP = q(g, "Sc4_X5_d2", "BART-CP", "bias"),
  gauss_Sc5d2_bias_LRC = q(g, "Sc5_d2", "LRC-BART-100", "bias"), gauss_Sc5d2_rmse_LRC = q(g, "Sc5_d2", "LRC-BART-100", "rmse"),
  gauss_Sc5d2_map = q(g, "Sc5_d2", "LRC-BART-100", "map_Rc"), gauss_Sc5d2_ess_real = q(g, "Sc5_d2", "LRC-BART-100", "ess_realized"),
  gauss_Sc5d2_bias_CBART = q(g, "Sc5_d2", "C-BART", "bias"), gauss_Sc5d2_bias_BARTCP = q(g, "Sc5_d2", "BART-CP", "bias"),
  gauss_Sc5d1_map = q(g, "Sc5_d1", "LRC-BART-100", "map_Rc"),
  surv_Sc1_rmse_LRC = q(s, "Sc1", "LRC-BART-100", "rmse", "med"), surv_Sc1_rmse_PP = q(s, "Sc1", "BART-PP", "rmse", "med"),
  surv_Sc1_bias_NP = q(s, "Sc1", "BART-NP", "bias", "med"), surv_Sc1_rmse_NP = q(s, "Sc1", "BART-NP", "rmse", "med"),
  surv_Sc1_bias_NP_rmst = q(s, "Sc1", "BART-NP", "bias", "rmst"),
  surv_Sc3_bias_LRC = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(s, sc, "LRC-BART-100", "bias", "med")),
  surv_Sc3_bias_PP = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(s, sc, "BART-PP", "bias", "med")),
  surv_Sc3_bias_NP = sapply(c("Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5"), function(sc) q(s, sc, "BART-NP", "bias", "med")),
  surv_Sc3_rmse_LRC = q(s, "Sc3_rho0", "LRC-BART-100", "rmse", "med"), surv_Sc3_rmse_PP = q(s, "Sc3_rho0", "BART-PP", "rmse", "med"),
  surv_Sc3_ctrl_bias_LRC = q(s, "Sc3_rho0", "LRC-BART-100", "ctrl_bias", "med"), surv_Sc3_ctrl_bias_PP = q(s, "Sc3_rho0", "BART-PP", "ctrl_bias", "med"),
  surv_Sc3_rmst_bias_LRC = q(s, "Sc3_rho0", "LRC-BART-100", "bias", "rmst"), surv_Sc3_rmst_bias_PP = q(s, "Sc3_rho0", "BART-PP", "bias", "rmst"),
  surv_Sc3_rmst_rmse_LRC = q(s, "Sc3_rho0", "LRC-BART-100", "rmse", "rmst"), surv_Sc3_rmst_rmse_PP = q(s, "Sc3_rho0", "BART-PP", "rmse", "rmst"),
  surv_Sc3_bias_CBART = q(s, "Sc3_rho0", "C-BART", "bias", "med"), surv_Sc3_bias_BARTCP = q(s, "Sc3_rho0", "BART-CP", "bias", "med"),
  surv_Sc4d1_map = c(q(s, "Sc4_X5_d1", "LRC-BART-100", "map_R", "med"), q(s, "Sc4_X5_d1", "LRC-BART-100", "map_Rc", "med")),
  surv_Sc4d1_bias_LRC = q(s, "Sc4_X5_d1", "LRC-BART-100", "bias", "med"), surv_Sc4d1_bias_CBART = q(s, "Sc4_X5_d1", "C-BART", "bias", "med"),
  surv_Sc5_bias_LRC = q(s, "Sc5_d1", "LRC-BART-100", "bias", "med"), surv_Sc5_map = q(s, "Sc5_d1", "LRC-BART-100", "map_Rc", "med"),
  surv_Sc5_bias_CBART = q(s, "Sc5_d1", "C-BART", "bias", "med"), surv_Sc5_bias_BARTCP = q(s, "Sc5_d1", "BART-CP", "bias", "med"),
  surv_Sc1_pow_LRC = q(s, "Sc1", "LRC-BART-100", "pow", "med"), surv_Sc1_pow_PP = q(s, "Sc1", "BART-PP", "pow", "med"),
  surv_Sc2_ceiling = q(s, "Sc2", "LRC-BART-100", "ess_ceiling", "med"), surv_Sc1_ceiling = q(s, "Sc1", "LRC-BART-100", "ess_ceiling", "med")
)
sink(file.path(OUT, "quoted_numbers.txt"))
print(quotes)
cat("\nCAHB: 02-validation/cahb/summary_table.csv\n"); print(cahb, row.names=FALSE)
cat("\nRegional ESS: 02-validation/sc4_boundary/regional_ess_summary.csv\n"); print(regional, row.names=FALSE)
cat("\nControl-profile surface: 02-validation/sc4_boundary/summary.csv\n"); print(boundary, row.names=FALSE)
cat("\nPaired control-profile contrasts: 02-validation/sc4_boundary/paired.csv\n"); print(boundary_paired, row.names=FALSE)
sink()
cat("tables written to", OUT, "\n")

# Unnumbered supplementary panels preserve the established S1-S6 numbering.
ca <- cahb[cahb$method %in% c("CAHB", "CAHB-p10", "CAHB-lit"), ]
ca_rows <- sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
  scen_short[ca$scenario], method_name[ca$method], f3(ca$bias), f3(ca$sd), f3(ca$rmse), f2(ca$map_R), f2(ca$map_Rc))
writeLines(c("\\begin{center}\\scriptsize", "\\begin{tabular}{llccccc}", "\\toprule",
  "Scenario & Transfer & Bias & SD & RMSE & Share, $\\mathcal R$ & Share, $\\mathcal R^c$ \\\\",
  "\\midrule", ca_rows, "\\bottomrule", "\\end{tabular}", "\\end{center}"), file.path(OUT, "tab_cahb_transfer.tex"))
bp <- boundary_paired[boundary_paired$metric %in% c("rmse_R", "rmse_Rc", "rmse_R_far", "rmse_Rc_far"), ]
region_name <- c(rmse_R="$X_5>2$", rmse_Rc="$X_5\\le2$", rmse_R_far="$X_5>2.5$", rmse_Rc_far="$X_5<1.5$")
br <- sprintf("%s & %s & %s & %s & %s & %s \\\\", scen_short[bp$scenario], region_name[bp$metric],
  f3(bp$lrc), f3(bp$pp), f3(bp$diff), f3(bp$se_paired))
writeLines(c("\\begin{center}\\scriptsize", "\\begin{tabular}{llcccc}", "\\toprule",
  "Scenario & Control profiles & LRC-BART & BART-PP & Difference & Paired SE \\\\",
  "\\midrule", br, "\\bottomrule", "\\end{tabular}", "\\end{center}"), file.path(OUT, "tab_boundary_detail.tex"))
