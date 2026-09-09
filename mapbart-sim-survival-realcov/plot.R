rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
mainDir <- "/Users/oliviazhang/Desktop/"

p_obs <- 10
size_option <- "default"   # "default" (RCT trt 200) or "small" (RCT trt 30)
size_suffix <- if (size_option == "small") "_n30" else "_n200"
# HierAFT (parametric) keeps a fixed prior; MAP-BART (mBART) is now
# calibrated per N_target multiplier and result files use _N<target>_
# suffixes (mirrors mapbart-sim-gaussian/plot.R).
aftv4_prior_vals  <- c(0.05, 0.5)
# AFTv2 power-prior RWD weights to display (files tagged _w<rwd_w>).
# Mirrors the rwd_w sweep in run_all.R: w=1 plus matched-N weights derived from
# mbart_target_Ns / n_rwd_nominal.  run_all.R overrides this at runtime.
aftv2_rwd_w_vals  <- if (size_option == "small") c(1, 0.3, 0.2, 0.1, 0.0733, 0.05) else c(1)
n_rwd_nominal     <- 300   # nominal RWD control count; labels shown as N = round(rwd_w * n_rwd_nominal)
# MAP-BART targets are the N_target multipliers (1.00, 0.75, 0.50) of the size's
# calibration N_target, plus the near-zero discrepancy reference (s^2=0.0001,
# maximal borrowing) tagged N=Inf.  default: N_target=200 -> 200/150/100;
# small: N_target=30 -> 30/22/15.  Listed in decreasing N so the MAP-BART
# curves/legend order as N=Inf, then descending.
mbart_target_Ns  <- if (size_option == "small") c(Inf, 90, 60, 30, 22, 15) else c(Inf, 200, 150, 100)
# MAP-BART control-arm BART configs to PLOT (must match configs mBART.R ran).
# run_all.R overrides this; ALL ntree:k pairs are overlaid as separate
# MAP-BART(nt,k,N=) curves.  With a single pair the label drops the nt,k prefix.
mbart_control_configs <- "5:2,10:2,50:2"
.mb_pairs <- strsplit(mbart_control_configs, "[, ]+")[[1]]
.mb_cfgs  <- lapply(.mb_pairs, function(p) {
  v <- strsplit(p, ":")[[1]]; nt <- as.integer(v[1]); kk <- as.numeric(v[2])
  list(nt = nt, k = kk,
       tag = paste0("_nt", nt, "_k", kk),                              # filename tag
       lab = if (length(.mb_pairs) > 1) sprintf("nt%d,k%g,", nt, kk) else "")  # method-label prefix
})
# Backward-compatibility alias used by HierAFT loops below.
prior_vals       <- aftv4_prior_vals

# Helper function to add separator rows between groups
add_subsc_separators <- function(df, subsc_col = "SubScenario_Label") {
  if (nrow(df) == 0 || !subsc_col %in% colnames(df)) return(df)

  # Get unique values in order of appearance
  unique_vals <- unique(as.character(df[[subsc_col]]))

  result <- data.frame()
  for (i in seq_along(unique_vals)) {
    group_data <- df[as.character(df[[subsc_col]]) == unique_vals[i], ]
    result <- rbind(result, group_data)

    # Add separator row after each group (except the last)
    if (i < length(unique_vals)) {
      sep_row <- group_data[1, ]
      sep_row[] <- ""
      result <- rbind(result, sep_row)
    }
  }

  return(result)
}

# =============================================================================
# build_rmst_column(sc)  --  RMST(tau)-ratio column (right side of the unified
# figure). Contains the SAME elements as the former standalone plot_rmst.R:
# RMST performance facet + summary table, arm-specific population RMST facet +
# tables, and (if present) subject-wise RMST facet + table. Reads the RMST
# columns fresh from res/*.RData. Returns a patchwork stack, or NULL if absent.
# =============================================================================
build_rmst_column <- function(sc, subsc = "E") {
  RMST_COLS <- c("c","iteration","rmst_tau","rmst_true","rmst_hat",
                 "bias_rmst","sd_rmst","rmse_rmst","w1distance_rmst","w2distance_rmst",
                 "ci_rmst","coverage_rmst","tp_rmst","tp_calibrated_rmst",
                 "bias.trt.rmst.pop","sd.trt.rmst.pop","w2distance.trt.rmst.pop",
                 "bias.ctrl.rmst.pop","sd.ctrl.rmst.pop","w2distance.ctrl.rmst.pop",
                 "bias_subj_rmst","pehe_subj_rmst")
  resdir <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/")
  read_rmst <- function(file, method) {
    if (!file.exists(file)) return(NULL)
    d <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.null(d) || !"bias_rmst" %in% names(d)) return(NULL)
    out <- d[, intersect(RMST_COLS, names(d)), drop = FALSE]; out$Method <- method; out
  }
  build_one <- function(sc) {
    suf <- function(extra = "") paste0("_sc", sc, subsc, if (sc == 3) "_cor0.7" else "", extra, "_alternative.RData")
    fr <- list()
    for (.aw in aftv2_rwd_w_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"AFTv2_p",p_obs,size_suffix,
                                        sub("\\.RData$", paste0("_w",.aw,".RData"), suf())),
                                 paste0("AFTv2(N=", round(.aw * n_rwd_nominal), ")"))))
    fr <- c(fr, list(read_rmst(paste0(resdir,"BARTv2_p",p_obs,size_suffix,suf()),"BARTv2")))
    for (pr in aftv4_prior_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"HierAFT_p",p_obs,size_suffix,suf(paste0("_prior",pr))), paste0("HierAFT(",pr,")"))))
    for (cf in .mb_cfgs) for (N in mbart_target_Ns)
      fr <- c(fr, list(read_rmst(paste0(resdir,"MAP-BART_p",p_obs,size_suffix,cf$tag,suf(paste0("_N",N))), paste0("MAP-BART(",cf$lab,"N=",N,")"))))
    fr <- fr[!sapply(fr, is.null)]; if (!length(fr)) return(NULL)
    out <- bind_rows(fr); out
  }
  dat <- build_one(sc)
  if (is.null(dat) || !nrow(dat)) return(NULL)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  grp <- "Method"
  # Main table: SAME columns as the median main summary table. RMSE = sqrt(mean(.))
  # like the median table; metrics fall back to "-" if absent in old result files.
  for (.cc in c("rmse_rmst","w1distance_rmst","w2distance_rmst","tp_calibrated_rmst",
                "sd.trt.rmst.pop","sd.ctrl.rmst.pop"))
    if (!.cc %in% names(dat)) dat[[.cc]] <- NA_real_
  .fmt  <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", m) }
  .fmtR <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", sqrt(m)) }
  rmst_main <- dat %>% group_by(Method) %>%
    summarise(Bias=.fmt(bias_rmst), SD=.fmt(sd_rmst), RMSE=.fmtR(rmse_rmst),
              W1Distance=.fmt(w1distance_rmst), W2Distance=.fmt(w2distance_rmst),
              CI_length=.fmt(ci_rmst), CI_coverage=.fmt(coverage_rmst),
              Power=.fmt(tp_rmst), Power_calib=.fmt(tp_calibrated_rmst),
              N=n(), N_bias_ge1=sum(abs(bias_rmst)>=1,na.rm=TRUE), .groups="drop") %>%
    mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)), .before=1)
  table_grob <- gridExtra::tableGrob(add_subsc_separators(rmst_main, "SubScenario_Label"), rows = NULL)
  long <- dat %>% pivot_longer(c(bias_rmst,sd_rmst,ci_rmst,coverage_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
  long$metric <- factor(long$metric, levels=c("bias_rmst","sd_rmst","ci_rmst","coverage_rmst"), labels=c("Bias","SD","CI length","Coverage"))
  thm <- theme_minimal() + theme(plot.title=element_text(size=14,face="bold",hjust=0.5), axis.text.x=element_text(size=9,angle=45,hjust=1), legend.position="none", strip.text=element_text(size=11,face="bold"), strip.background=element_rect(fill="gray90",color="gray50"), panel.border=element_rect(color="black",fill=NA,linewidth=0.5))
  p <- ggplot(long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("RMST(τ)-ratio performance — Scenario ",sc," (τ = ",sprintf("%.1f",mean(dat$rmst_tau,na.rm=TRUE)),")"), x=NULL, y=NULL) + thm
  p <- p + facet_wrap(~metric,scales="free_y",nrow=1)
  arm_long <- dat %>% pivot_longer(c(bias.trt.rmst.pop,w2distance.trt.rmst.pop,bias.ctrl.rmst.pop,w2distance.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control")), Metric=ifelse(grepl("^bias",metric),"Bias","W2Distance"))
  # Arm tables: SAME columns as the median population-median tables (Group, SubScenario_Label, Method, Bias, W2Distance).
  arm_summary <- arm_long %>% group_by(Method, Arm, Metric) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=Metric, values_from=m)
  # SD column: mean posterior SD of the per-arm estimate (matches the ATE table's
  # SD = mean(sd_rmst); uses per-replicate posterior SD cols sd.{trt,ctrl}.rmst.pop).
  arm_sd <- dat %>% pivot_longer(c(sd.trt.rmst.pop,sd.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control"))) %>%
    group_by(Method, Arm) %>% summarise(SD=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop")
  arm_summary <- dplyr::left_join(arm_summary, arm_sd, by=c("Method","Arm"))
  arm_tab_trt <- arm_summary %>% filter(Arm=="Treatment") %>% dplyr::select(-Arm) %>% mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)), Group="Treatment") %>% dplyr::select(Group, SubScenario_Label, Method, Bias, SD, W2Distance)
  arm_tab_ctrl <- arm_summary %>% filter(Arm=="Control") %>% dplyr::select(-Arm) %>% mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)), Group="Control") %>% dplyr::select(Group, SubScenario_Label, Method, Bias, SD, W2Distance)
  table_grob_arm <- wrap_plots(gridExtra::tableGrob(add_subsc_separators(arm_tab_trt,"SubScenario_Label"),rows=NULL), gridExtra::tableGrob(add_subsc_separators(arm_tab_ctrl,"SubScenario_Label"),rows=NULL), ncol=2)
  p_arm <- ggplot(arm_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("Arm-specific population RMST(τ) — Scenario ",sc,"  (estimate vs DGP truth, per arm)"), x=NULL,y=NULL) + thm
  p_arm <- p_arm + facet_grid(Arm~Metric,scales="free_y")
  # Table-slot height grows with the number of methods so all rows show in full.
  .tbl_h <- max(2, length(unique(as.character(dat$Method))) * 0.16)
  has_subj <- all(c("bias_subj_rmst","pehe_subj_rmst") %in% names(dat))
  if (has_subj) {
    subj_long <- dat %>% pivot_longer(c(bias_subj_rmst,pehe_subj_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
    subj_long$metric <- factor(subj_long$metric, levels=c("bias_subj_rmst","pehe_subj_rmst"), labels=c("Bias","PEHE"))
    # Subject-wise table: SAME columns as the median subject-wise table (SubScenario_Label, Method, Bias, PEHE).
    subj_tab <- subj_long %>% group_by(Method, metric) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=metric, values_from=m) %>%
      mutate(SubScenario_Label=factor(paste0("SubSc ",subsc))) %>% dplyr::select(SubScenario_Label, Method, Bias, PEHE)
    table_grob_subj <- gridExtra::tableGrob(add_subsc_separators(subj_tab,"SubScenario_Label"), rows = NULL)
    p_subj <- ggplot(subj_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
      geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
      labs(title=paste0("Subject-wise RMST(τ)-ratio — Scenario ",sc,"  (per-subject estimate vs truth)"), x=NULL,y=NULL) + thm
    p_subj <- p_subj + facet_wrap(~metric,scales="free_y",nrow=1)
    # Blank bottom section (matches the median column's variance/sigma section) so
    # the two columns have identical row layout and align horizontally.
    return(p / table_grob / p_arm / table_grob_arm / p_subj / table_grob_subj /
             patchwork::plot_spacer() / patchwork::plot_spacer() +
             plot_layout(heights = c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h, 4, .tbl_h)))
  }
  p / table_grob / p_arm / table_grob_arm /
    patchwork::plot_spacer() / patchwork::plot_spacer() +
    plot_layout(heights = c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h))
}

for (sc in c(1, 2, 3)) {
tryCatch({

if (sc == 1 | sc == 2 | sc == 3){

  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # Loop through sub-scenarios (E only)
  subscenarios <- c("E")

  # sc3 uses unmeasured-confounding data tagged _cor<rho> (rho = 0.7); sc1/2 have no cor tag.
  .cor_sfx <- if (sc == 3) "_cor0.7" else ""

  for (subsc in subscenarios) {

    # Initialize empty lists for this sub-scenario
    subsc_res_list_ATE <- list()
    subsc_res_list_sigma <- list()

    # Construct file suffix based on subsc
    if (subsc == "NULL") {
      file_suffix <- paste0("_sc", sc, .cor_sfx, "_alternative.RData")
      subsc_label <- "NULL"
    } else {
      file_suffix <- paste0("_sc", sc, subsc, "", .cor_sfx, "_alternative.RData")
      subsc_label <- subsc
    }

    # Try to read AFTv1 results
    tryCatch({
      AFTv1_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv1_p",p_obs,size_suffix,file_suffix))
      AFTv1_res_ATE <- AFTv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv1_res_sigma <- AFTv1_res[, c("c","iteration",
                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                "pehe_subj","bias_subj")]
      AFTv1_res_ATE$Method <- "AFTv1"
      AFTv1_res_sigma$Method <- "AFTv1"
      AFTv1_res_ATE$SubScenario <- subsc_label
      AFTv1_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv1_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv1_res_sigma
    }, error = function(e) {})

    # Try to read AFTv2 results for each rwd_w value (files tagged _w<rwd_w>)
    for (rwd_w_val in aftv2_rwd_w_vals) {
      tryCatch({
        aftv2_file <- sub("\\.RData$", paste0("_w", rwd_w_val, ".RData"), paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv2_p",p_obs,size_suffix,file_suffix))
        if (file.exists(aftv2_file)) {
          AFTv2_res <- readRDS(aftv2_file)
          AFTv2_res_ATE <- AFTv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv2_res_sigma <- AFTv2_res[, c("c","iteration",
                                      "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                      "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                      "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                      "pehe_subj","bias_subj")]
          aftv2_label <- paste0("AFTv2(N=", round(rwd_w_val * n_rwd_nominal), ")")
          AFTv2_res_ATE$Method <- aftv2_label
          AFTv2_res_sigma$Method <- aftv2_label
          AFTv2_res_ATE$SubScenario <- subsc_label
          AFTv2_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv2_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv2_res_sigma
        }
      }, error = function(e) {})
    }

    # Try to read AFTv3 results
    tryCatch({
      AFTv3_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv3_p",p_obs,size_suffix,file_suffix))
      AFTv3_res_ATE <- AFTv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv3_res_sigma <- AFTv3_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      AFTv3_res_ATE$Method <- "AFTv3"
      AFTv3_res_sigma$Method <- "AFTv3"
      AFTv3_res_ATE$SubScenario <- subsc_label
      AFTv3_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv3_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv3_res_sigma
    }, error = function(e) {})

    # Try to read HierAFT results for each prior value
    for (prior_val in prior_vals) {
      tryCatch({
        aftv4_file_suffix <- paste0("_sc", sc, subsc, "", .cor_sfx, "_prior", prior_val, "_alternative.RData")
        file_path <- paste0(mainDir,"/mapbart-sim-survival-realcov/res/HierAFT_p",p_obs,size_suffix,aftv4_file_suffix)
        if (file.exists(file_path)) {
          HierAFT_res <- readRDS(file_path)
          HierAFT_res_ATE <- HierAFT_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          HierAFT_res_sigma <- HierAFT_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          HierAFT_res_ATE$Method <- paste0("HierAFT(", prior_val, ")")
          HierAFT_res_sigma$Method <- paste0("HierAFT(", prior_val, ")")
          HierAFT_res_ATE$SubScenario <- subsc_label
          HierAFT_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- HierAFT_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- HierAFT_res_sigma
        }
      }, error = function(e) {})
    }

    # Try to read BARTv1 / BARTv2 results
    for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
    tryCatch({
      BART_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival-realcov/res/",bver,"_p",p_obs,size_suffix,file_suffix))
      BART_res <- BART_res[BART_res$alpha == 0.95 & BART_res$beta == 2, ]
      BART_res_ATE <- BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      BART_res_sigma <- BART_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      BART_res_ATE$Method <- bver
      BART_res_sigma$Method <- bver
      BART_res_ATE$SubScenario <- subsc_label
      BART_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- BART_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- BART_res_sigma
    }, error = function(e) {})
    }

    # Try to read MAP-BART results for each control-config x N_target value
    for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
      tryCatch({
        mbart_file_suffix <- paste0(cf$tag, "_sc", sc, subsc, "", .cor_sfx, "_N", target_N, "_alternative.RData")
        file_path <- paste0(mainDir,"/mapbart-sim-survival-realcov/res/MAP-BART_p",p_obs,size_suffix,mbart_file_suffix)
        if (file.exists(file_path)) {
          MAP_BART_res <- readRDS(file_path)
          MAP_BART_res <- MAP_BART_res[MAP_BART_res$alpha == 0.95 & MAP_BART_res$beta == 2, ]
          MAP_BART_res_ATE <- MAP_BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          MAP_BART_res_sigma <- MAP_BART_res[, c("c","iteration",
                                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma","bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                "bias.trt.median.pop","sd.trt.median.pop","w2distance.trt.median.pop",
                                                "bias.ctrl.median.pop","sd.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                "pehe_subj","bias_subj")]
          MAP_BART_res_ATE$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
          MAP_BART_res_sigma$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
          MAP_BART_res_ATE$SubScenario <- subsc_label
          MAP_BART_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- MAP_BART_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- MAP_BART_res_sigma
        }
      }, error = function(e) {})
    }

    # Combine results for this sub-scenario (only if we have data)
    if (length(subsc_res_list_ATE) > 0) {
      subsc_res_ATE <- do.call(rbind, subsc_res_list_ATE)
      all_res_ATE <- rbind(all_res_ATE, subsc_res_ATE)
    }

    if (length(subsc_res_list_sigma) > 0) {
      subsc_res_sigma <- do.call(rbind, subsc_res_list_sigma)
      all_res_sigma <- rbind(all_res_sigma, subsc_res_sigma)
    }
  }

  # Assign to res and res_sigma for plotting
  res <- all_res_ATE
  res_sigma <- all_res_sigma

  # Load Type I error from null hypothesis files (if current hypo is not null)
  all_null_FP <- data.frame()
  for (subsc in subscenarios) {
    # Construct null file suffix
    if (subsc == "NULL") {
      null_file_suffix <- paste0("_sc", sc, .cor_sfx, "_null.RData")
      subsc_label <- "NULL"
    } else {
      null_file_suffix <- paste0("_sc", sc, subsc, "", .cor_sfx, "_null.RData")
      subsc_label <- subsc
    }

    # Try to load null results for each method
    methods_list <- list(
      list(name = "AFTv1", file_prefix = "AFTv1"),
      list(name = "AFTv2", file_prefix = "AFTv2"),
      list(name = "AFTv3", file_prefix = "AFTv3"),
      list(name = "BARTv1", file_prefix = "BARTv1"),
      list(name = "BARTv2", file_prefix = "BARTv2"),
      list(name = "BARTv3", file_prefix = "BARTv3")
    )

    for (method_info in methods_list) {
      tryCatch({
        null_file <- paste0(mainDir,"/mapbart-sim-survival-realcov/res/", method_info$file_prefix, "_p", p_obs, size_suffix, null_file_suffix)
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)

          # For BART, filter by alpha and beta
          if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3")) {
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          }

          # Calculate mean FP (Type I error)
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              SubScenario = subsc_label,
              Method = method_info$name,
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {
        # Silently skip if null file doesn't exist
      })
    }

    # Try to load HierAFT null results for each prior value
    for (prior_val in prior_vals) {
      tryCatch({
        aftv4_null_suffix <- paste0("_sc", sc, subsc, .cor_sfx, "_prior", prior_val, "_null.RData")
        null_file <- paste0(mainDir,"/mapbart-sim-survival-realcov/res/HierAFT_p", p_obs, size_suffix, aftv4_null_suffix)
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              SubScenario = subsc_label,
              Method = paste0("HierAFT(", prior_val, ")"),
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {})
    }

    # Try to load MAP-BART null results for each control-config x N_target value
    for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
      tryCatch({
        mbart_null_suffix <- paste0(cf$tag, "_sc", sc, subsc, "", .cor_sfx, "_N", target_N, "_null.RData")
        null_file <- paste0(mainDir,"/mapbart-sim-survival-realcov/res/MAP-BART_p", p_obs, size_suffix, mbart_null_suffix)
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)
          null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              SubScenario = subsc_label,
              Method = paste0("MAP-BART(", cf$lab, "N=", target_N, ")"),
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {})
    }
  }
}

# Check if we have any data to plot
if (nrow(res) == 0) {
  stop("No results data found. Please check that result files exist and paths are correct.")
}

# Create SubScenario labels for scenarios 1, 2, and 3
if (sc == 1 | sc == 2 | sc == 3) {
  if (nrow(res) > 0 && "SubScenario" %in% colnames(res)) {
    res$SubScenario_Label <- factor(paste0("SubSc ", res$SubScenario),
                                     levels = c("SubSc E"))
  }
  if (nrow(res_sigma) > 0 && "SubScenario" %in% colnames(res_sigma)) {
    res_sigma$SubScenario_Label <- factor(paste0("SubSc ", res_sigma$SubScenario),
                                           levels = c("SubSc E"))
  }
}

# Create dynamic method levels including AFTv2, HierAFT and MAP-BART with N labels.
aftv2_methods <- paste0("AFTv2(N=", round(aftv2_rwd_w_vals * n_rwd_nominal), ")")
aftv4_methods <- paste0("HierAFT(",     aftv4_prior_vals, ")")
mbart_methods <- unlist(lapply(.mb_cfgs, function(cf) paste0("MAP-BART(", cf$lab, "N=", mbart_target_Ns, ")")))
all_method_levels <- c("AFTv1", aftv2_methods, "AFTv3", "BARTv1", "BARTv2", "BARTv3", aftv4_methods, mbart_methods)

if (nrow(res) > 0) {
  res$Method <- factor(res$Method, levels = all_method_levels)
}
if (nrow(res_sigma) > 0) {
  res_sigma$Method <- factor(res_sigma$Method, levels = all_method_levels)
}

# Define specific colors for each method (default ggplot2 colors)
n_methods <- length(all_method_levels)
default_colors <- scales::hue_pal()(n_methods)
method_colors <- setNames(default_colors, all_method_levels)

if (sc == 1 | sc == 2 | sc == 3){

  summary_table <- res %>%
    filter(!is.na(bias)) %>%
    group_by(SubScenario_Label, SubScenario, Method) %>%
    summarise(
      Bias = sprintf("%.2f", mean(bias, na.rm = TRUE)),
      SD = sprintf("%.2f", mean(sd, na.rm = TRUE)),
      RMSE = sprintf("%.2f", sqrt(mean(rmse, na.rm = TRUE))),
      W1Distance = sprintf("%.2f", mean(w1distance, na.rm = TRUE)),
      W2Distance = sprintf("%.2f", mean(w2distance, na.rm = TRUE)),
      CI_length = sprintf("%.2f", mean(ci, na.rm = TRUE)),
      CI_coverage = sprintf("%.2f", mean(coverage, na.rm = TRUE)),
      Power = sprintf("%.2f", mean(tp, na.rm = TRUE)),
      Power_calib = sprintf("%.2f", mean(tp_calibrated, na.rm = TRUE)),
      N = n(),
      N_bias_ge1 = sum(abs(bias) >= 1, na.rm = TRUE),
      .groups = "drop"
    )

  # Add Type I error from null files if available
  if (nrow(all_null_FP) > 0) {
    summary_table <- summary_table %>%
      left_join(all_null_FP, by = c("SubScenario", "Method")) %>%
      rename(`Type I error` = Type_I_error)
  }

  # Remove SubScenario column (keep only SubScenario_Label)
  summary_table <- summary_table %>% dplyr::select(-SubScenario)

  # Add separator rows between subscenarios and create table
  summary_table_sep <- add_subsc_separators(summary_table, "SubScenario_Label")
  table_grob <- gridExtra::tableGrob(summary_table_sep, rows = NULL)

  res <- res %>%
    filter(!is.na(bias)) %>%
    pivot_longer(cols = c(bias, sd, rmse, w1distance, w2distance),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))

  facet_plot <- ggplot(res, aes(x = Method, y = value, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Treatment Effect Estimation Performance Across Sub-Scenarios",
      x = NULL,
      y = NULL
    ) +
    facet_grid(SubScenario_Label ~ metric, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray50"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

  # Variance estimation plots (only if we have sigma data)
  if (nrow(res_sigma) > 0) {
    res_sigma_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
      pivot_longer(cols = c(bias.trt.sigma, sd.trt.sigma, w2distance.trt.sigma,
                            bias.ctrl.sigma, sd.ctrl.sigma, w2distance.ctrl.sigma),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.sigma", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "sd" ~ "SD",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_sigma <- res_sigma_long %>%
      group_by(SubScenario_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      # Reorder columns: Bias, SD, W2Distance
      dplyr::select(SubScenario_Label, Method, group, Bias, SD, W2Distance)

    # Split into treatment and control tables
    summary_table_sigma_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_sigma_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_sigma_trt_sep <- add_subsc_separators(summary_table_sigma_trt, "SubScenario_Label")
    summary_table_sigma_ctrl_sep <- add_subsc_separators(summary_table_sigma_ctrl, "SubScenario_Label")
    table_grob_sigma_trt <- gridExtra::tableGrob(summary_table_sigma_trt_sep, rows = NULL)
    table_grob_sigma_ctrl <- gridExtra::tableGrob(summary_table_sigma_ctrl_sep, rows = NULL)
    table_grob_sigma <- wrap_plots(table_grob_sigma_trt, table_grob_sigma_ctrl, ncol = 2)

    facet_plot_sigma <- ggplot(res_sigma_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Variance Estimation Performance Across Sub-Scenarios",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + SubScenario_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # Check if subject-wise columns exist
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
      # Subject-wise bias and RMSE plots
      res_subj_long <- res_sigma %>%
        filter(Method %in% all_method_levels) %>%
        pivot_longer(cols = c(bias_subj, pehe_subj),
                     names_to = "metric",
                     values_to = "value") %>%
        mutate(
          metric_type = case_when(
            metric == "bias_subj" ~ "Bias",
            metric == "pehe_subj" ~ "PEHE",
            TRUE ~ metric
          )
        )

      if (nrow(res_subj_long) > 0) {
        summary_table_subj <- res_subj_long %>%
          group_by(SubScenario_Label, Method, metric_type) %>%
          summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
          pivot_wider(names_from = metric_type, values_from = Mean)

        # Add separator rows between subscenarios and create table
        summary_table_subj_sep <- add_subsc_separators(summary_table_subj, "SubScenario_Label")
        table_grob_subj <- gridExtra::tableGrob(summary_table_subj_sep, rows = NULL)

        facet_plot_subj <- ggplot(res_subj_long, aes(x = Method, y = value, fill = Method)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
          geom_boxplot(alpha = 0.7) +
          stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
          scale_fill_manual(values = method_colors) +
          labs(
            title = "Subject-wise Performance Across Sub-Scenarios",
            x = NULL,
            y = NULL
          ) +
          facet_grid(SubScenario_Label ~ metric_type, scales = "free_y") +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray90", color = "gray50"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
          )
      } else {
        facet_plot_subj <- NULL
        table_grob_subj <- NULL
      }
    } else {
      facet_plot_subj <- NULL
      table_grob_subj <- NULL
    }

    # Population median plots
    res_pop_median_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
      pivot_longer(cols = c(bias.trt.median.pop, w2distance.trt.median.pop,
                            bias.ctrl.median.pop, w2distance.ctrl.median.pop),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.median\\.pop", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_pop_median <- res_pop_median_long %>%
      group_by(SubScenario_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean)

    # SD column: mean posterior SD of the per-arm estimate (matches the ATE table's
    # SD = mean(sd); uses the per-replicate posterior SD cols sd.{trt,ctrl}.median.pop).
    sd_pop_median <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
      pivot_longer(cols = c(sd.trt.median.pop, sd.ctrl.median.pop),
                   names_to = "metric", values_to = "value") %>%
      mutate(group = factor(ifelse(grepl("trt", metric), "Treatment", "Control"),
                            levels = c("Treatment", "Control"))) %>%
      group_by(SubScenario_Label, Method, group) %>%
      summarise(SD = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop")
    summary_table_pop_median <- dplyr::left_join(summary_table_pop_median, sd_pop_median,
      by = c("SubScenario_Label", "Method", "group"))

    # Split into treatment and control tables (Bias, SD, W2Distance)
    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      mutate(Group = "Treatment") %>%
      dplyr::select(Group, SubScenario_Label, Method, Bias, SD, W2Distance)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      mutate(Group = "Control") %>%
      dplyr::select(Group, SubScenario_Label, Method, Bias, SD, W2Distance)

    # Add separator rows between subscenarios and create tables
    summary_table_pop_median_trt_sep <- add_subsc_separators(summary_table_pop_median_trt, "SubScenario_Label")
    summary_table_pop_median_ctrl_sep <- add_subsc_separators(summary_table_pop_median_ctrl, "SubScenario_Label")
    table_grob_pop_median_trt <- gridExtra::tableGrob(summary_table_pop_median_trt_sep, rows = NULL)
    table_grob_pop_median_ctrl <- gridExtra::tableGrob(summary_table_pop_median_ctrl_sep, rows = NULL)
    table_grob_pop_median <- wrap_plots(table_grob_pop_median_trt, table_grob_pop_median_ctrl, ncol = 2)

    facet_plot_pop_median <- ggplot(res_pop_median_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Median Performance",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + SubScenario_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # Table-slot height grows with the number of methods so all rows show in full.
    # NOTE: `res` is long-format here (melted into `metric`/`value`), so count
    # distinct methods directly rather than via a (now-absent) `bias` column.
    .n_meth <- length(unique(as.character(res$Method)))
    .tbl_h  <- max(2, .n_meth * 0.16)
    # Combine ATE, population median, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      .heights <- c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h, 4, .tbl_h)
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = .heights)
    } else {
      .heights <- c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h)
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = .heights)
    }
  } else {
    # Only ATE plot if no sigma data
    .heights <- c(4, .tbl_h <- max(2, length(unique(as.character(res$Method))) * 0.16))
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = .heights)
  }

  # ---- Append the RMST(tau)-ratio column (right) beside the median-ratio column (left) ----
  rmst_col <- build_rmst_column(sc)
  if (!is.null(rmst_col)) {
    final_plot <- final_plot | rmst_col
  }

  # Dynamically calculate dimensions based on number of sub-scenarios
  n_subsc <- length(unique(res$SubScenario_Label))
  n_plots <- if (nrow(res_sigma) > 0) 4 else 1  # ATE, subj, pop_median, sigma or just ATE

  # Calculate total dimensions.  Height tracks the layout's relative heights
  # (~3 in per unit) so the taller, multi-config tables render in full.
  total_height <- if (exists(".heights")) sum(.heights) * n_subsc * 3 + 10 else n_subsc * 15 * n_plots + 10
  total_width <- if (!is.null(rmst_col)) 36 else 18  # double width: median (left) + RMST (right)

  ggsave(paste0(file.path(mainDir),"/mapbart-sim-survival-realcov/inserts/p",p_obs,size_suffix,"_sc",sc,"_all_subsc_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}

}, error = function(e) {
  message(sprintf("Skipped sc = %d: %s", sc, conditionMessage(e)))
})
}
