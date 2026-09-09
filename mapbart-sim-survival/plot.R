rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
mainDir <- "/Users/oliviazhang/Desktop/"

p_obs <- 10
# HierAFT (parametric) keeps a fixed prior; MAP-BART (mBART) is now
# calibrated per N_target multiplier and result files use _N<target>_
# suffixes (mirrors mapbart-sim-gaussian/plot.R).
aftv4_prior_vals <- 0.25
# Standalone default; run_all.R overrides this from the ess calibration file.
# Listed in decreasing N (Inf first) for MAP-BART curve/legend order.
mbart_target_Ns  <- c(Inf, 100, 75, 50)
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
# figure). RMST performance facet + summary table, arm-specific population RMST
# facet + tables, and (if present) subject-wise RMST facet + table. Reads the
# RMST columns fresh from res/*.RData. Returns a patchwork stack, or NULL if
# absent. Uses `subsc` from the enclosing scope (matches the median column).
# =============================================================================
build_rmst_column <- function(sc) {
  RMST_COLS <- c("c","iteration","rmst_tau","rmst_true","rmst_hat",
                 "bias_rmst","sd_rmst","rmse_rmst","w1distance_rmst","w2distance_rmst",
                 "ci_rmst","coverage_rmst","tp_rmst","tp_calibrated_rmst",
                 "bias.trt.rmst.pop","w2distance.trt.rmst.pop",
                 "bias.ctrl.rmst.pop","w2distance.ctrl.rmst.pop",
                 "bias_subj_rmst","pehe_subj_rmst")
  resdir <- paste0(mainDir, "/mapbart-sim-survival/res/")
  read_rmst <- function(file, method) {
    if (!file.exists(file)) return(NULL)
    d <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.null(d) || !"bias_rmst" %in% names(d)) return(NULL)
    out <- d[, intersect(RMST_COLS, names(d)), drop = FALSE]; out$Method <- method; out
  }
  build_one <- function(sc, cor_tag = NULL) {
    cortxt <- if (is.null(cor_tag)) "" else paste0("_cor", cor_tag)
    suf <- function(extra = "") paste0("_sc", sc, subsc, cortxt, extra, "_alternative.RData")
    fr <- list(read_rmst(paste0(resdir,"AFTv1_p",p_obs,suf()),"AFTv1"),
               read_rmst(paste0(resdir,"AFTv2_p",p_obs,suf()),"AFTv2"),
               read_rmst(paste0(resdir,"AFTv3_p",p_obs,suf()),"AFTv3"),
               read_rmst(paste0(resdir,"BARTv1_p",p_obs,suf()),"BARTv1"),
               read_rmst(paste0(resdir,"BARTv2_p",p_obs,suf()),"BARTv2"),
               read_rmst(paste0(resdir,"BARTv3_p",p_obs,suf()),"BARTv3"))
    for (pr in aftv4_prior_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"HierAFT_p",p_obs,suf(paste0("_prior",pr))), paste0("HierAFT(",pr,")"))))
    for (N in mbart_target_Ns)
      fr <- c(fr, list(read_rmst(paste0(resdir,"MAP-BART_p",p_obs,suf(paste0("_N",N))), paste0("MAP-BART(N=",N,")"))))
    fr <- fr[!sapply(fr, is.null)]; if (!length(fr)) return(NULL)
    out <- bind_rows(fr); out$cor <- if (is.null(cor_tag)) NA_character_ else as.character(cor_tag); out
  }
  dat <- if (sc %in% c(1,2)) build_one(sc) else bind_rows(lapply(c(-0.5,0,0.5), function(cv) build_one(3,cv)))
  if (is.null(dat) || !nrow(dat)) return(NULL)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  grp <- if (sc == 3) c("Method","cor") else "Method"
  # Finalize a summary table the SAME way as the median-ratio tables: insert
  # separator rows, and for sc3 keep a plain "SubSc E" label plus a separate
  # Correlation column, separating rows by the SubScenario+rho combination
  # (the combined "..., rho=" label is used only to place separators, then dropped).
  fin_tab <- function(tb, label_col = "SubScenario_Label") {
    if (sc == 3) {
      tb <- tb %>% arrange(.data[[label_col]], Correlation, Method) %>%
        mutate(.sep = paste0(.data[[label_col]], ", rho=", Correlation))
      add_subsc_separators(tb, ".sep") %>% dplyr::select(-.sep)
    } else add_subsc_separators(tb, label_col)
  }
  for (.cc in c("rmse_rmst","w1distance_rmst","w2distance_rmst","tp_calibrated_rmst"))
    if (!.cc %in% names(dat)) dat[[.cc]] <- NA_real_
  .fmt  <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", m) }
  .fmtR <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", sqrt(m)) }
  rmst_main <- dat %>% group_by(across(all_of(grp))) %>%
    summarise(Bias=.fmt(bias_rmst), SD=.fmt(sd_rmst), RMSE=.fmtR(rmse_rmst),
              W1Distance=.fmt(w1distance_rmst), W2Distance=.fmt(w2distance_rmst),
              CI_length=.fmt(ci_rmst), CI_coverage=.fmt(coverage_rmst),
              Power=.fmt(tp_rmst), Power_calib=.fmt(tp_calibrated_rmst),
              N=n(), N_bias_ge1=sum(abs(bias_rmst)>=1,na.rm=TRUE), .groups="drop") %>%
    mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)), .before=1)
  if (sc==3) rmst_main <- rmst_main %>% rename(Correlation=cor) %>% relocate(Correlation, .after=SubScenario_Label)
  table_grob <- gridExtra::tableGrob(fin_tab(rmst_main), rows = NULL)
  long <- dat %>% pivot_longer(c(bias_rmst,sd_rmst,ci_rmst,coverage_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
  long$metric <- factor(long$metric, levels=c("bias_rmst","sd_rmst","ci_rmst","coverage_rmst"), labels=c("Bias","SD","CI length","Coverage"))
  thm <- theme_minimal() + theme(plot.title=element_text(size=14,face="bold",hjust=0.5), axis.text.x=element_text(size=9,angle=45,hjust=1), legend.position="none", strip.text=element_text(size=11,face="bold"), strip.background=element_rect(fill="gray90",color="gray50"), panel.border=element_rect(color="black",fill=NA,linewidth=0.5))
  p <- ggplot(long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("RMST(τ)-ratio performance — Scenario ",sc," (τ = ",sprintf("%.1f",mean(dat$rmst_tau,na.rm=TRUE)),")"), x=NULL, y=NULL) + thm
  p <- if (sc==3) p + facet_grid(cor~metric,scales="free_y") else p + facet_wrap(~metric,scales="free_y",nrow=1)
  arm_long <- dat %>% pivot_longer(c(bias.trt.rmst.pop,w2distance.trt.rmst.pop,bias.ctrl.rmst.pop,w2distance.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control")), Metric=ifelse(grepl("^bias",metric),"Bias","W2Distance"))
  arm_summary <- arm_long %>% group_by(across(all_of(c(grp, "Arm", "Metric")))) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=Metric, values_from=m)
  mk_arm <- function(g) {
    t <- arm_summary %>% filter(Arm==g) %>% dplyr::select(-Arm) %>% mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)), Group=g)
    if (sc==3) t %>% rename(Correlation=cor) %>% dplyr::select(Group, SubScenario_Label, Correlation, Method, Bias, W2Distance)
    else       t %>% dplyr::select(Group, SubScenario_Label, Method, Bias, W2Distance)
  }
  arm_tab_trt <- mk_arm("Treatment"); arm_tab_ctrl <- mk_arm("Control")
  table_grob_arm <- wrap_plots(gridExtra::tableGrob(fin_tab(arm_tab_trt),rows=NULL), gridExtra::tableGrob(fin_tab(arm_tab_ctrl),rows=NULL), ncol=2)
  p_arm <- ggplot(arm_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("Arm-specific population RMST(τ) — Scenario ",sc,"  (estimate vs DGP truth, per arm)"), x=NULL,y=NULL) + thm
  p_arm <- if (sc==3) p_arm + facet_grid(Arm+cor~Metric,scales="free_y") else p_arm + facet_grid(Arm~Metric,scales="free_y")
  has_subj <- all(c("bias_subj_rmst","pehe_subj_rmst") %in% names(dat))
  if (has_subj) {
    subj_long <- dat %>% pivot_longer(c(bias_subj_rmst,pehe_subj_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
    subj_long$metric <- factor(subj_long$metric, levels=c("bias_subj_rmst","pehe_subj_rmst"), labels=c("Bias","PEHE"))
    subj_tab <- subj_long %>% group_by(across(all_of(c(grp, "metric")))) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=metric, values_from=m) %>%
      mutate(SubScenario_Label=factor(paste0("SubSc ",subsc)))
    subj_tab <- if (sc==3) subj_tab %>% rename(Correlation=cor) %>% dplyr::select(SubScenario_Label, Correlation, Method, Bias, PEHE)
                else        subj_tab %>% dplyr::select(SubScenario_Label, Method, Bias, PEHE)
    table_grob_subj <- gridExtra::tableGrob(fin_tab(subj_tab), rows = NULL)
    p_subj <- ggplot(subj_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
      geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
      labs(title=paste0("Subject-wise RMST(τ)-ratio — Scenario ",sc,"  (per-subject estimate vs truth)"), x=NULL,y=NULL) + thm
    p_subj <- if (sc==3) p_subj + facet_grid(cor~metric,scales="free_y") else p_subj + facet_wrap(~metric,scales="free_y",nrow=1)
    return(p / table_grob / p_arm / table_grob_arm / p_subj / table_grob_subj /
             patchwork::plot_spacer() / patchwork::plot_spacer() +
             plot_layout(heights = c(4, 2, 4, 2, 4, 1.5, 4, 2)))
  }
  p / table_grob / p_arm / table_grob_arm /
    patchwork::plot_spacer() / patchwork::plot_spacer() +
    plot_layout(heights = c(4, 2, 4, 2, 4, 2))
}

for (sc in 1:3) {
tryCatch({

if (sc == 1 | sc == 2){

  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # Loop through sub-scenarios (E only)
  subscenarios <- c("E")

  for (subsc in subscenarios) {

    # Initialize empty lists for this sub-scenario
    subsc_res_list_ATE <- list()
    subsc_res_list_sigma <- list()

    # Construct file suffix based on subsc
    if (subsc == "NULL") {
      file_suffix <- paste0("_sc", sc, "_alternative.RData")
      subsc_label <- "NULL"
    } else {
      file_suffix <- paste0("_sc", sc, subsc, "_alternative.RData")
      subsc_label <- subsc
    }

    # Try to read AFTv1 results
    tryCatch({
      AFTv1_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival/res/AFTv1_p",p_obs,file_suffix))
      AFTv1_res_ATE <- AFTv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv1_res_sigma <- AFTv1_res[, c("c","iteration",
                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                "pehe_subj","bias_subj")]
      AFTv1_res_ATE$Method <- "AFTv1"
      AFTv1_res_sigma$Method <- "AFTv1"
      AFTv1_res_ATE$SubScenario <- subsc_label
      AFTv1_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv1_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv1_res_sigma
    }, error = function(e) {})

    # Try to read AFTv2 results
    tryCatch({
      AFTv2_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival/res/AFTv2_p",p_obs,file_suffix))
      AFTv2_res_ATE <- AFTv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv2_res_sigma <- AFTv2_res[, c("c","iteration",
                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                  "pehe_subj","bias_subj")]
      AFTv2_res_ATE$Method <- "AFTv2"
      AFTv2_res_sigma$Method <- "AFTv2"
      AFTv2_res_ATE$SubScenario <- subsc_label
      AFTv2_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv2_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv2_res_sigma
    }, error = function(e) {})

    # Try to read AFTv3 results
    tryCatch({
      AFTv3_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival/res/AFTv3_p",p_obs,file_suffix))
      AFTv3_res_ATE <- AFTv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv3_res_sigma <- AFTv3_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
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
        aftv4_file_suffix <- paste0("_sc", sc, subsc, "_prior", prior_val, "_alternative.RData")
        file_path <- paste0(mainDir,"/mapbart-sim-survival/res/HierAFT_p",p_obs,aftv4_file_suffix)
        if (file.exists(file_path)) {
          HierAFT_res <- readRDS(file_path)
          HierAFT_res_ATE <- HierAFT_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          HierAFT_res_sigma <- HierAFT_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
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
      BART_res <- readRDS(paste0(mainDir,"/mapbart-sim-survival/res/",bver,"_p",p_obs,file_suffix))
      BART_res <- BART_res[BART_res$alpha == 0.95 & BART_res$beta == 2, ]
      BART_res_ATE <- BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      BART_res_sigma <- BART_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      BART_res_ATE$Method <- bver
      BART_res_sigma$Method <- bver
      BART_res_ATE$SubScenario <- subsc_label
      BART_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- BART_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- BART_res_sigma
    }, error = function(e) {})
    }

    # Try to read MAP-BART results for each N_target value
    for (target_N in mbart_target_Ns) {
      tryCatch({
        mbart_file_suffix <- paste0("_sc", sc, subsc, "_N", target_N, "_alternative.RData")
        file_path <- paste0(mainDir,"/mapbart-sim-survival/res/MAP-BART_p",p_obs,mbart_file_suffix)
        if (file.exists(file_path)) {
          MAP_BART_res <- readRDS(file_path)
          MAP_BART_res <- MAP_BART_res[MAP_BART_res$alpha == 0.95 & MAP_BART_res$beta == 2, ]
          MAP_BART_res_ATE <- MAP_BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          MAP_BART_res_sigma <- MAP_BART_res[, c("c","iteration",
                                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                "pehe_subj","bias_subj")]
          MAP_BART_res_ATE$Method <- paste0("MAP-BART(N=", target_N, ")")
          MAP_BART_res_sigma$Method <- paste0("MAP-BART(N=", target_N, ")")
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
      null_file_suffix <- paste0("_sc", sc, "_null.RData")
      subsc_label <- "NULL"
    } else {
      null_file_suffix <- paste0("_sc", sc, subsc, "_null.RData")
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
        null_file <- paste0(mainDir,"/mapbart-sim-survival/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
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
        aftv4_null_suffix <- paste0("_sc", sc, subsc, "_prior", prior_val, "_null.RData")
        null_file <- paste0(mainDir,"/mapbart-sim-survival/res/HierAFT_p", p_obs, aftv4_null_suffix)
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

    # Try to load MAP-BART null results for each N_target value
    for (target_N in mbart_target_Ns) {
      tryCatch({
        mbart_null_suffix <- paste0("_sc", sc, subsc, "_N", target_N, "_null.RData")
        null_file <- paste0(mainDir,"/mapbart-sim-survival/res/MAP-BART_p", p_obs, mbart_null_suffix)
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)
          null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              SubScenario = subsc_label,
              Method = paste0("MAP-BART(N=", target_N, ")"),
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {})
    }
  }
}
if (sc == 3){
  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # Loop through sub-scenarios (E only) and correlations
  subscenarios <- c("E")
  cor_vals <- c(-0.5, 0, 0.5)

  # Track which files were found
  files_checked <- list()
  files_found <- list()

  for (subsc in subscenarios) {
    for (cor_val in cor_vals) {

      # Initialize empty lists for this sub-scenario and correlation
      subsc_res_list_ATE <- list()
      subsc_res_list_sigma <- list()

      # Construct file suffix based on subsc
      # For sc == 3: include "cor" prefix before correlation value
      if (subsc == "NULL") {
        file_suffix <- paste0("_sc", sc, "_cor", cor_val, "_alternative.RData")
        subsc_label <- "NULL"
      } else {
        file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_alternative.RData")
        subsc_label <- subsc
      }

      # Try to read AFTv1 results
      tryCatch({
        file_path <- paste0(mainDir,"mapbart-sim-survival/res/AFTv1_p",p_obs,file_suffix)
        files_checked[[length(files_checked) + 1]] <- file_path
        if (file.exists(file_path)) {
          AFTv1_res <- readRDS(file_path)
          files_found[[length(files_found) + 1]] <- file_path
          AFTv1_res_ATE <- AFTv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv1_res_sigma <- AFTv1_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv1_res_ATE$Method <- "AFTv1"
          AFTv1_res_sigma$Method <- "AFTv1"
          AFTv1_res_ATE$SubScenario <- subsc_label
          AFTv1_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv1_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv1_res_sigma
        }
      }, error = function(e) {})

      # Try to read AFTv2 results
      tryCatch({
        file_path <- paste0(mainDir,"mapbart-sim-survival/res/AFTv2_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          AFTv2_res <- readRDS(file_path)
          AFTv2_res_ATE <- AFTv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv2_res_sigma <- AFTv2_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv2_res_ATE$Method <- "AFTv2"
          AFTv2_res_sigma$Method <- "AFTv2"
          AFTv2_res_ATE$SubScenario <- subsc_label
          AFTv2_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv2_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv2_res_sigma
        }
      }, error = function(e) {})

      # Try to read AFTv3 results
      tryCatch({
        file_path <- paste0(mainDir,"mapbart-sim-survival/res/AFTv3_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          AFTv3_res <- readRDS(file_path)
          AFTv3_res_ATE <- AFTv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv3_res_sigma <- AFTv3_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv3_res_ATE$Method <- "AFTv3"
          AFTv3_res_sigma$Method <- "AFTv3"
          AFTv3_res_ATE$SubScenario <- subsc_label
          AFTv3_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- AFTv3_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- AFTv3_res_sigma
        }
      }, error = function(e) {})

      # Try to read HierAFT results for each prior value
      for (prior_val in prior_vals) {
        tryCatch({
          aftv4_file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_prior", prior_val, "_alternative.RData")
          file_path <- paste0(mainDir,"mapbart-sim-survival/res/HierAFT_p",p_obs,aftv4_file_suffix)
          if (file.exists(file_path)) {
            HierAFT_res <- readRDS(file_path)
            HierAFT_res_ATE <- HierAFT_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            HierAFT_res_sigma <- HierAFT_res[, c("c","iteration",
                                          "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                          "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                          "bias.trt.median.pop","w2distance.trt.median.pop",
                                          "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
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
        file_path <- paste0(mainDir,"mapbart-sim-survival/res/",bver,"_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          BART_res <- readRDS(file_path)
          BART_res <- BART_res[BART_res$alpha == 0.95 & BART_res$beta == 2, ]
          BART_res_ATE <- BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          BART_res_sigma <- BART_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          BART_res_ATE$Method <- bver
          BART_res_sigma$Method <- bver
          BART_res_ATE$SubScenario <- subsc_label
          BART_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- BART_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- BART_res_sigma
        }
      }, error = function(e) {})
      }

      # Try to read MAP-BART results for each N_target value
      for (target_N in mbart_target_Ns) {
        tryCatch({
          mbart_file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_N", target_N, "_alternative.RData")
          file_path <- paste0(mainDir,"mapbart-sim-survival/res/MAP-BART_p",p_obs,mbart_file_suffix)
          if (file.exists(file_path)) {
            MAP_BART_res <- readRDS(file_path)
            MAP_BART_res <- MAP_BART_res[MAP_BART_res$alpha == 0.95 & MAP_BART_res$beta == 2, ]
            MAP_BART_res_ATE <- MAP_BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            MAP_BART_res_sigma <- MAP_BART_res[, c("c","iteration",
                                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                  "pehe_subj","bias_subj")]
            MAP_BART_res_ATE$Method <- paste0("MAP-BART(N=", target_N, ")")
            MAP_BART_res_sigma$Method <- paste0("MAP-BART(N=", target_N, ")")
            MAP_BART_res_ATE$SubScenario <- subsc_label
            MAP_BART_res_sigma$SubScenario <- subsc_label
            subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- MAP_BART_res_ATE
            subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- MAP_BART_res_sigma
          }
        }, error = function(e) {})
      }

      # Combine results for this sub-scenario and correlation (only if we have data)
      if (length(subsc_res_list_ATE) > 0) {
        subsc_res_ATE <- do.call(rbind, subsc_res_list_ATE)
        all_res_ATE <- rbind(all_res_ATE, subsc_res_ATE)
      }

      if (length(subsc_res_list_sigma) > 0) {
        subsc_res_sigma <- do.call(rbind, subsc_res_list_sigma)
        all_res_sigma <- rbind(all_res_sigma, subsc_res_sigma)
      }
    }
  }

  # Assign to res and res_sigma for plotting
  if (nrow(all_res_ATE) > 0) {
    res <- all_res_ATE
  } else {
    stop("No results files found for scenario 3")
  }

  if (nrow(all_res_sigma) > 0) {
    res_sigma <- all_res_sigma
  } else {
    res_sigma <- data.frame()
  }

  # Load Type I error from null hypothesis files
  all_null_FP <- data.frame()
  methods_list <- list(
    list(name = "AFTv1", file_prefix = "AFTv1"),
    list(name = "AFTv2", file_prefix = "AFTv2"),
    list(name = "AFTv3", file_prefix = "AFTv3"),
    list(name = "BARTv1", file_prefix = "BARTv1"),
    list(name = "BARTv2", file_prefix = "BARTv2"),
    list(name = "BARTv3", file_prefix = "BARTv3")
  )

  for (subsc in subscenarios) {
    for (cor_val in cor_vals) {
      # Construct null file suffix based on subsc
      if (subsc == "NULL") {
        null_file_suffix <- paste0("_sc", sc, "_cor", cor_val, "_null.RData")
        subsc_label <- "NULL"
      } else {
        null_file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_null.RData")
        subsc_label <- subsc
      }

      for (method_info in methods_list) {
        tryCatch({
          null_file <- paste0(mainDir,"mapbart-sim-survival/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
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
                c = cor_val,
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
          aftv4_null_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_prior", prior_val, "_null.RData")
          null_file <- paste0(mainDir,"mapbart-sim-survival/res/HierAFT_p", p_obs, aftv4_null_suffix)
          if (file.exists(null_file)) {
            null_res <- readRDS(null_file)
            if ("fp" %in% colnames(null_res)) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                SubScenario = subsc_label,
                c = cor_val,
                Method = paste0("HierAFT(", prior_val, ")"),
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {})
      }

      # Try to load MAP-BART null results for each N_target value
      for (target_N in mbart_target_Ns) {
        tryCatch({
          mbart_null_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_N", target_N, "_null.RData")
          null_file <- paste0(mainDir,"mapbart-sim-survival/res/MAP-BART_p", p_obs, mbart_null_suffix)
          if (file.exists(null_file)) {
            null_res <- readRDS(null_file)
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            if ("fp" %in% colnames(null_res)) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                SubScenario = subsc_label,
                c = cor_val,
                Method = paste0("MAP-BART(N=", target_N, ")"),
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {})
      }
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

# Create dynamic method levels including HierAFT and MAP-BART with prior values
# HierAFT still labels by prior value; MAP-BART labels by N_target.
aftv4_methods <- paste0("HierAFT(",     aftv4_prior_vals, ")")
mbart_methods <- paste0("MAP-BART(N=", mbart_target_Ns,  ")")
all_method_levels <- c("AFTv1", "AFTv2", "AFTv3", "BARTv1", "BARTv2", "BARTv3", aftv4_methods, mbart_methods)

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

if (sc == 1 | sc ==2){

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
        # Order group factor with Treatment first (on top)
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_sigma <- res_sigma_long %>%
      group_by(SubScenario_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      # Reorder columns: Bias, SD, W2Distance
      dplyr::select(SubScenario_Label, Method, group, Bias, SD, W2Distance)

    # Split into treatment and control tables with clear labels
    summary_table_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_trt_sep <- add_subsc_separators(summary_table_trt, "SubScenario_Label")
    summary_table_ctrl_sep <- add_subsc_separators(summary_table_ctrl, "SubScenario_Label")
    table_grob_trt <- gridExtra::tableGrob(summary_table_trt_sep, rows = NULL)
    table_grob_ctrl <- gridExtra::tableGrob(summary_table_ctrl_sep, rows = NULL)

    # Combine treatment and control tables horizontally
    table_grob_sigma <- wrap_plots(table_grob_trt, table_grob_ctrl, ncol = 2)

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

    # Split into treatment and control tables
    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

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

    # Combine ATE, population median, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 1.5, 4, 2))
    } else {
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 2))
    }
  } else {
    # Only ATE plot if no sigma data
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = c(4, 1))
  }

  # ---- Append the RMST(tau)-ratio column (right) beside the median-ratio column (left) ----
  rmst_col <- build_rmst_column(sc)
  if (!is.null(rmst_col)) {
    final_plot <- final_plot | rmst_col
  }

  # Dynamically calculate dimensions based on number of sub-scenarios
  n_subsc <- length(unique(res$SubScenario_Label))
  n_plots <- if (nrow(res_sigma) > 0) 4 else 1  # ATE, subj, pop_median, sigma or just ATE

  # Calculate total dimensions
  total_height <- n_subsc * 15 * n_plots + 10
  total_width <- if (!is.null(rmst_col)) 36 else 18  # double width: median (left) + RMST (right)

  ggsave(paste0(file.path(mainDir),"/mapbart-sim-survival/inserts/p",p_obs,"_sc",sc,"_all_subsc_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}
if (sc == 3){

  # Check if res has data and required columns
  if (nrow(res) == 0) {
    stop("No data remaining for sc == 3 after filtering. All iterations may have |bias| >= 1. Check the filtering threshold or data quality.")
  }
  if (!"SubScenario_Label" %in% colnames(res)) {
    stop("SubScenario_Label column not found in res. Check data loading.")
  }

  summary_table <- res %>%
    filter(!is.na(bias)) %>%
    group_by(SubScenario_Label, SubScenario, c, Method) %>%
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
    ) %>%
    # Order by subscenario first, then correlation
    arrange(SubScenario_Label, c, Method)

  # Add Type I error from null files if available
  if (nrow(all_null_FP) > 0) {
    summary_table <- summary_table %>%
      left_join(all_null_FP, by = c("SubScenario", "c", "Method")) %>%
      rename(`Type I error` = Type_I_error)
  }

  # Create combined SubScenario + Correlation label for separator function
  summary_table$SubSc_Cor_Label <- paste0(summary_table$SubScenario_Label, ", rho=", summary_table$c)

  # Remove SubScenario column (keep only SubScenario_Label) and rename c to Correlation
  summary_table <- summary_table %>% dplyr::select(-SubScenario) %>% rename(Correlation = c)

  # Add separator rows between each subscenario-correlation combination and create table
  summary_table_sep <- add_subsc_separators(summary_table, "SubSc_Cor_Label")
  # Remove the SubSc_Cor_Label column from display (it's redundant with SubScenario_Label + Correlation)
  summary_table_sep <- summary_table_sep %>% dplyr::select(-SubSc_Cor_Label)
  table_grob <- gridExtra::tableGrob(summary_table_sep, rows = NULL)

  # Create combined SubScenario + Correlation label for faceting
  # Order: group by subscenario first, then by correlation within each subscenario
  if (nrow(res) > 0) {
    res$SubSc_Cor_Label <- factor(paste0(res$SubScenario_Label, ", rho=", res$c),
                                   levels = as.vector(t(outer(c("SubSc E"),
                                                              c(-0.5, 0, 0.5),
                                                              function(s, r) paste0(s, ", rho=", r)))))
  }
  if (nrow(res_sigma) > 0) {
    res_sigma$SubSc_Cor_Label <- factor(paste0(res_sigma$SubScenario_Label, ", rho=", res_sigma$c),
                                         levels = as.vector(t(outer(c("SubSc E"),
                                                                    c(-0.5, 0, 0.5),
                                                                    function(s, r) paste0(s, ", rho=", r)))))
  }

  res <- res %>%
    filter(!is.na(bias)) %>%
    pivot_longer(cols = c(bias, sd, rmse, w1distance, w2distance),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))

  # Determine number of unique combinations that exist
  n_subsc_cor <- length(unique(res$SubSc_Cor_Label))

  facet_plot <- ggplot(res, aes(x = Method, y = value, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Treatment Effect Estimation Performance Across Sub-Scenarios and Correlations",
      x = NULL,
      y = NULL
    ) +
    facet_grid(SubSc_Cor_Label ~ metric, scales = "free_y") +
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
      group_by(SubSc_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      dplyr::select(SubSc_Cor_Label, c, Method, group, Bias, SD, W2Distance) %>%
      rename(Correlation = c) %>%
      # Order by subscenario first, then correlation
      arrange(SubSc_Cor_Label, Method)

    summary_table_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_trt_sep <- add_subsc_separators(summary_table_trt, "SubSc_Cor_Label")
    summary_table_ctrl_sep <- add_subsc_separators(summary_table_ctrl, "SubSc_Cor_Label")
    table_grob_trt <- gridExtra::tableGrob(summary_table_trt_sep, rows = NULL)
    table_grob_ctrl <- gridExtra::tableGrob(summary_table_ctrl_sep, rows = NULL)
    table_grob_sigma <- wrap_plots(table_grob_trt, table_grob_ctrl, ncol = 2)

    facet_plot_sigma <- ggplot(res_sigma_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Variance Estimation Performance Across Sub-Scenarios and Correlations",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + SubSc_Cor_Label ~ metric_type, scales = "free_y") +
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

    # Subject-wise plots
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
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
          group_by(SubSc_Cor_Label, c, Method, metric_type) %>%
          summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
          pivot_wider(names_from = metric_type, values_from = Mean) %>%
          rename(Correlation = c) %>%
          # Order by subscenario first, then correlation
          arrange(SubSc_Cor_Label, Method)

        # Add separator rows between subscenarios and create table
        summary_table_subj_sep <- add_subsc_separators(summary_table_subj, "SubSc_Cor_Label")
        table_grob_subj <- gridExtra::tableGrob(summary_table_subj_sep, rows = NULL)

        facet_plot_subj <- ggplot(res_subj_long, aes(x = Method, y = value, fill = Method)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
          geom_boxplot(alpha = 0.7) +
          stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
          scale_fill_manual(values = method_colors) +
          labs(
            title = "Subject-wise Performance Across Sub-Scenarios and Correlations",
            x = NULL,
            y = NULL
          ) +
          facet_grid(SubSc_Cor_Label ~ metric_type, scales = "free_y") +
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
      group_by(SubSc_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      rename(Correlation = c) %>%
      # Order by subscenario first, then correlation
      arrange(SubSc_Cor_Label, Method)

    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_pop_median_trt_sep <- add_subsc_separators(summary_table_pop_median_trt, "SubSc_Cor_Label")
    summary_table_pop_median_ctrl_sep <- add_subsc_separators(summary_table_pop_median_ctrl, "SubSc_Cor_Label")
    table_grob_pop_median_trt <- gridExtra::tableGrob(summary_table_pop_median_trt_sep, rows = NULL)
    table_grob_pop_median_ctrl <- gridExtra::tableGrob(summary_table_pop_median_ctrl_sep, rows = NULL)
    table_grob_pop_median <- wrap_plots(table_grob_pop_median_trt, table_grob_pop_median_ctrl, ncol = 2)

    facet_plot_pop_median <- ggplot(res_pop_median_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Median Performance Across Sub-Scenarios and Correlations",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + SubSc_Cor_Label ~ metric_type, scales = "free_y") +
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

    # Combine ATE, population median, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 1.5, 4, 2))
    } else {
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 2))
    }
  } else {
    # Only ATE plot if no sigma data
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = c(4, 1))
  }

  # ---- Append the RMST(tau)-ratio column (right) beside the median-ratio column (left) ----
  rmst_col <- build_rmst_column(sc)
  if (!is.null(rmst_col)) {
    final_plot <- final_plot | rmst_col
  }

  # Dynamically calculate dimensions based on number of sub-scenarios and correlations
  n_subsc <- length(unique(res$SubSc_Cor_Label))
  n_plots <- if (nrow(res_sigma) > 0) 4 else 1  # ATE, subj, pop_median, sigma or just ATE

  # Calculate total dimensions
  total_height <- n_subsc * 12 * n_plots + 10
  total_width <- if (!is.null(rmst_col)) 36 else 18  # double width: median (left) + RMST (right)

  ggsave(paste0(file.path(mainDir),"mapbart-sim-survival/inserts/p",p_obs,"_sc",sc,"_all_subsc_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}

}, error = function(e) {
  message(sprintf("Skipped sc = %d: %s", sc, conditionMessage(e)))
})
}
