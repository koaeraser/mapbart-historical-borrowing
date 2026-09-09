rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
mainDir <- "/Users/oliviazhang/Desktop/"

p_obs <- 10

# Single-arm toggle (set by run_all.R). When FALSE, the saved plot filenames get
# a "_singlearm" tag so single-arm results don't overwrite the full-control ones.
RCT_ctrl <- FALSE
rct_ctrl_tag <- ""

# HierLM discrepancy-variance prior sweep.
lmv4_prior_vals <- c(0.05, 0.5)

# LMv2 power-prior RWD weights: derived below from mbart_target_Ns so LMv2 shows
# every N MAP-BART targets (rwd_w = N / n_rwd_nominal) plus 1 (full borrowing).
n_rwd_nominal   <- 300   # nominal RWD control count; labels shown as N = round(rwd_w * n_rwd_nominal)

# MAP and mBART now key on N_target rather than the raw s^2 value: each
# analysis file loops over multiple N_target multipliers from its
# calibration .RData and writes one result file per target (with the
# integer target N in the filename).  List the N values to read here.
map_target_Ns   <- c(200, 150, 100)
# Calibrated s^2_star targets are the N_target multipliers (1.00, 0.75, 0.50)
# of the calibration's N_target=200, plus the near-zero discrepancy reference
# (s^2=0.0001, maximal borrowing) tagged N=Inf.  Listed in decreasing N so the
# MAP-BART curves/legend order as N=Inf, 200, 150, 100.
mbart_target_Ns <- c(Inf, 200, 150, 100)
# LMv2 runs at every N MAP-BART targets (rwd_w = N/n_rwd_nominal) plus full borrowing (rwd_w=1).
lmv2_rwd_w_vals <- sort(unique(c(1, round(mbart_target_Ns[is.finite(mbart_target_Ns)] / n_rwd_nominal, 4))), decreasing = TRUE)
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
      file_suffix <- paste0("_sc", sc, "_alternative", rct_ctrl_tag, ".RData")
      subsc_label <- "NULL"
    } else {
      file_suffix <- paste0("_sc", sc, subsc, "_alternative", rct_ctrl_tag, ".RData")
      subsc_label <- subsc
    }

    # Try to read LMv1 results
    tryCatch({
      LMv1_res <- readRDS(paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv1_p",p_obs,file_suffix))
      LMv1_res_ATE <- LMv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv1_res_sigma <- LMv1_res[, c("c","iteration",
                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                "pehe_subj","bias_subj")]
      LMv1_res_ATE$Method <- "LMv1"
      LMv1_res_sigma$Method <- "LMv1"
      LMv1_res_ATE$SubScenario <- subsc_label
      LMv1_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv1_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv1_res_sigma
    }, error = function(e) {
    })

    # Try to read LMv2 results for each rwd_w value (files tagged _w<rwd_w>)
    for (rwd_w_val in lmv2_rwd_w_vals) {
      tryCatch({
        if (subsc == "NULL") {
          lmv2_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv2_p",p_obs,"_sc",sc,"_alternative_w",rwd_w_val,rct_ctrl_tag,".RData")
        } else {
          lmv2_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv2_p",p_obs,"_sc",sc,subsc,"_alternative_w",rwd_w_val,rct_ctrl_tag,".RData")
        }
        if (file.exists(lmv2_file)) {
          LMv2_res <- readRDS(lmv2_file)
          LMv2_res_ATE <- LMv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv2_res_sigma <- LMv2_res[, c("c","iteration",
                                      "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                      "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                      "bias.trt.median.pop","w2distance.trt.median.pop",
                                      "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                      "pehe_subj","bias_subj")]
          lmv2_label <- paste0("LMv2(N=", round(rwd_w_val * n_rwd_nominal), ")")
          LMv2_res_ATE$Method <- lmv2_label
          LMv2_res_sigma$Method <- lmv2_label
          LMv2_res_ATE$SubScenario <- subsc_label
          LMv2_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv2_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv2_res_sigma
        }
      }, error = function(e) {
      })
    }

    # Try to read LMv3 results
    tryCatch({
      LMv3_res <- readRDS(paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv3_p",p_obs,file_suffix))
      LMv3_res_ATE <- LMv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv3_res_sigma <- LMv3_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      LMv3_res_ATE$Method <- "LMv3"
      LMv3_res_sigma$Method <- "LMv3"
      LMv3_res_ATE$SubScenario <- subsc_label
      LMv3_res_sigma$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv3_res_ATE
      subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv3_res_sigma
    }, error = function(e) {
    })

    # Try to read HierLM results for each prior value
    for (prior_val in lmv4_prior_vals) {
      tryCatch({
        # Construct file path with prior
        if (subsc == "NULL") {
          lmv4_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/HierLM_p",p_obs,"_sc",sc,"_prior",prior_val,"_alternative", rct_ctrl_tag, ".RData")
        } else {
          lmv4_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/HierLM_p",p_obs,"_sc",sc,subsc,"_prior",prior_val,"_alternative", rct_ctrl_tag, ".RData")
        }

        if (file.exists(lmv4_file)) {
          HierLM_res <- readRDS(lmv4_file)
          HierLM_res_ATE <- HierLM_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          HierLM_res_sigma <- HierLM_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          HierLM_res_ATE$Method <- paste0("HierLM(", prior_val, ")")
          HierLM_res_sigma$Method <- paste0("HierLM(", prior_val, ")")
          HierLM_res_ATE$SubScenario <- subsc_label
          HierLM_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- HierLM_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- HierLM_res_sigma
        }
      }, error = function(e) {
      })
    }

    # Try to read MAP results for each N_target value (ATE only, no sigma data)
    for (target_N in map_target_Ns) {
      tryCatch({
        # Construct file path with target N
        if (subsc == "NULL") {
          map_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP_p",p_obs,"_sc",sc,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
        } else {
          map_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP_p",p_obs,"_sc",sc,subsc,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
        }

        if (file.exists(map_file)) {
          MAP_res <- readRDS(map_file)
          MAP_res_ATE <- MAP_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          MAP_res_sigma <- MAP_res[, c("c","iteration",
                                      "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                      "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                      "bias.trt.median.pop","w2distance.trt.median.pop",
                                      "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                      "pehe_subj","bias_subj")]
          MAP_res_ATE$Method <- paste0("MAP(N=", target_N, ")")
          MAP_res_sigma$Method <- paste0("MAP(N=", target_N, ")")
          MAP_res_ATE$SubScenario <- subsc_label
          MAP_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- MAP_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- MAP_res_sigma
        }
      }, error = function(e) {
      })
    }

    # Try to read PSCL results (ATE only, no sigma data, may have NA values)
    # Keep all rows including NAs - plotting functions will handle with na.rm = TRUE
    tryCatch({
      PSCL_res <- readRDS(paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/PSCL_p",p_obs,file_suffix))
      PSCL_res_ATE <- PSCL_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      PSCL_res_ATE$Method <- "PSCL"
      PSCL_res_ATE$SubScenario <- subsc_label
      subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- PSCL_res_ATE
      # PSCL doesn't have sigma data, so no sigma entry
    }, error = function(e) {
    })

    # Try to read BARTv1 / BARTv2 results
    for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
    tryCatch({
      BART_res <- readRDS(paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/",bver,"_p",p_obs,file_suffix))
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
    }, error = function(e) {
    })
    }

    # Try to read mBART results for each control-config x N_target value
    for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
      tryCatch({
        # Construct file path with control-config tag and target N
        if (subsc == "NULL") {
          mbart_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,cf$tag,"_sc",sc,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
        } else {
          mbart_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,cf$tag,"_sc",sc,subsc,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
        }

        if (file.exists(mbart_file)) {
          mBART_res <- readRDS(mbart_file)
          mBART_res <- mBART_res[mBART_res$alpha == 0.95 & mBART_res$beta == 2, ]
          mBART_res_ATE <- mBART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          mBART_res_sigma <- mBART_res[, c("c","iteration",
                                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                "pehe_subj","bias_subj")]
          mBART_res_ATE$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
          mBART_res_sigma$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
          mBART_res_ATE$SubScenario <- subsc_label
          mBART_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- mBART_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- mBART_res_sigma
        }
      }, error = function(e) {
        # Silently skip if file doesn't exist
      })
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
      null_file_suffix <- paste0("_sc", sc, "_null", rct_ctrl_tag, ".RData")
      subsc_label <- "NULL"
    } else {
      null_file_suffix <- paste0("_sc", sc, subsc, "_null", rct_ctrl_tag, ".RData")
      subsc_label <- subsc
    }

    # Try to load null results for each method
    methods_list <- list(
      list(name = "LMv1", file_prefix = "LMv1"),
      list(name = "LMv2", file_prefix = "LMv2"),
      list(name = "LMv3", file_prefix = "LMv3"),
      list(name = "HierLM", file_prefix = "HierLM"),
      list(name = "MAP", file_prefix = "MAP"),
      list(name = "PSCL", file_prefix = "PSCL"),
      list(name = "BARTv1", file_prefix = "BARTv1"),
      list(name = "BARTv2", file_prefix = "BARTv2"),
      list(name = "BARTv3", file_prefix = "BARTv3")
    )

    for (method_info in methods_list) {
      tryCatch({
        null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)

          # For BART, filter by alpha and beta
          if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3")) {
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          }

          # Calculate mean FP (Type I error) - use na.rm = TRUE for PSCL which may have NAs
          if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
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

    # Try to load MAP-BART null results for each control-config x N_target value
    for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
      tryCatch({
        if (subsc == "NULL") {
          null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, cf$tag, "_sc", sc, "_N", target_N, "_null", rct_ctrl_tag, ".RData")
        } else {
          null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, cf$tag, "_sc", sc, subsc, "_N", target_N, "_null", rct_ctrl_tag, ".RData")
        }
        if (file.exists(null_file)) {
          null_res <- readRDS(null_file)
          null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              SubScenario = subsc_label,
              Method = paste0("MAP-BART(", cf$lab, "N=", target_N, ")"),
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {
        # Silently skip if null file doesn't exist
      })
    }
  }
}
if (sc == 3){
  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # Loop through sub-scenarios (E only) and correlations
  subscenarios <- c("E")
  cor_vals <- 0.7

  for (subsc in subscenarios) {
    for (cor_val in cor_vals) {

      # Initialize empty lists for this sub-scenario and correlation
      subsc_res_list_ATE <- list()
      subsc_res_list_sigma <- list()

      # Construct file suffix based on subsc
      # For sc == 3: include "cor" prefix before correlation value
      if (subsc == "NULL") {
        file_suffix <- paste0("_sc", sc, "_cor", cor_val, "_alternative", rct_ctrl_tag, ".RData")
        subsc_label <- "NULL"
      } else {
        file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_alternative", rct_ctrl_tag, ".RData")
        subsc_label <- subsc
      }

      # Try to read LMv1 results
      tryCatch({
        file_path <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv1_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          LMv1_res <- readRDS(file_path)
          LMv1_res_ATE <- LMv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv1_res_sigma <- LMv1_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          LMv1_res_ATE$Method <- "LMv1"
          LMv1_res_sigma$Method <- "LMv1"
          LMv1_res_ATE$SubScenario <- subsc_label
          LMv1_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv1_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv1_res_sigma
        }
      }, error = function(e) {
      })

      # Try to read LMv2 results for each rwd_w value (files tagged _w<rwd_w>)
      for (rwd_w_val in lmv2_rwd_w_vals) {
        tryCatch({
          if (subsc == "NULL") {
            lmv2_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv2_p",p_obs,"_sc",sc,"_cor",cor_val,"_alternative_w",rwd_w_val,rct_ctrl_tag,".RData")
          } else {
            lmv2_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv2_p",p_obs,"_sc",sc,subsc,"_cor",cor_val,"_alternative_w",rwd_w_val,rct_ctrl_tag,".RData")
          }
          if (file.exists(lmv2_file)) {
            LMv2_res <- readRDS(lmv2_file)
            LMv2_res_ATE <- LMv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            LMv2_res_sigma <- LMv2_res[, c("c","iteration",
                                          "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                          "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                          "bias.trt.median.pop","w2distance.trt.median.pop",
                                          "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                          "pehe_subj","bias_subj")]
            lmv2_label <- paste0("LMv2(N=", round(rwd_w_val * n_rwd_nominal), ")")
            LMv2_res_ATE$Method <- lmv2_label
            LMv2_res_sigma$Method <- lmv2_label
            LMv2_res_ATE$SubScenario <- subsc_label
            LMv2_res_sigma$SubScenario <- subsc_label
            subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv2_res_ATE
            subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv2_res_sigma
          }
        }, error = function(e) {
        })
      }

      # Try to read LMv3 results
      tryCatch({
        file_path <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/LMv3_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          LMv3_res <- readRDS(file_path)
          LMv3_res_ATE <- LMv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv3_res_sigma <- LMv3_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          LMv3_res_ATE$Method <- "LMv3"
          LMv3_res_sigma$Method <- "LMv3"
          LMv3_res_ATE$SubScenario <- subsc_label
          LMv3_res_sigma$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- LMv3_res_ATE
          subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- LMv3_res_sigma
        }
      }, error = function(e) {
      })

      # Try to read HierLM results for each prior value
      for (prior_val in lmv4_prior_vals) {
        tryCatch({
          # Construct file path with prior
          if (subsc == "NULL") {
            lmv4_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/HierLM_p",p_obs,"_sc",sc,"_cor",cor_val,"_prior",prior_val,"_alternative", rct_ctrl_tag, ".RData")
          } else {
            lmv4_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/HierLM_p",p_obs,"_sc",sc,subsc,"_cor",cor_val,"_prior",prior_val,"_alternative", rct_ctrl_tag, ".RData")
          }

          if (file.exists(lmv4_file)) {
            HierLM_res <- readRDS(lmv4_file)
            HierLM_res_ATE <- HierLM_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            HierLM_res_sigma <- HierLM_res[, c("c","iteration",
                                          "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                          "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                          "bias.trt.median.pop","w2distance.trt.median.pop",
                                          "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                          "pehe_subj","bias_subj")]
            HierLM_res_ATE$Method <- paste0("HierLM(", prior_val, ")")
            HierLM_res_sigma$Method <- paste0("HierLM(", prior_val, ")")
            HierLM_res_ATE$SubScenario <- subsc_label
            HierLM_res_sigma$SubScenario <- subsc_label
            subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- HierLM_res_ATE
            subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- HierLM_res_sigma
          }
        }, error = function(e) {
        })
      }

      # Try to read MAP results for each N_target value (ATE only, no sigma data)
      for (target_N in map_target_Ns) {
        tryCatch({
          # Construct file path with target N
          if (subsc == "NULL") {
            map_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP_p",p_obs,"_sc",sc,"_cor",cor_val,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
          } else {
            map_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP_p",p_obs,"_sc",sc,subsc,"_cor",cor_val,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
          }

          if (file.exists(map_file)) {
            MAP_res <- readRDS(map_file)
            MAP_res_ATE <- MAP_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            MAP_res_sigma <- MAP_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
            MAP_res_ATE$Method <- paste0("MAP(N=", target_N, ")")
            MAP_res_sigma$Method <- paste0("MAP(N=", target_N, ")")
            MAP_res_ATE$SubScenario <- subsc_label
            MAP_res_sigma$SubScenario <- subsc_label
            subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- MAP_res_ATE
            subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- MAP_res_sigma
          }
        }, error = function(e) {
        })
      }

      # Try to read PSCL results (ATE only, no sigma data, may have NA values)
      # Keep all rows including NAs - plotting functions will handle with na.rm = TRUE
      tryCatch({
        file_path <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/PSCL_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          PSCL_res <- readRDS(file_path)
          PSCL_res_ATE <- PSCL_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          PSCL_res_ATE$Method <- "PSCL"
          PSCL_res_ATE$SubScenario <- subsc_label
          subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- PSCL_res_ATE
          # PSCL doesn't have sigma data, so no sigma entry
        }
      }, error = function(e) {
      })

      # Try to read BARTv1 / BARTv2 results
      for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
      tryCatch({
        file_path <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/",bver,"_p",p_obs,file_suffix)
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
      }, error = function(e) {
      })
      }

      # Try to read mBART results for each control-config x N_target value
      for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
        tryCatch({
          # Construct file path with control-config tag and target N
          if (subsc == "NULL") {
            mbart_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,cf$tag,"_sc",sc,"_cor",cor_val,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
          } else {
            mbart_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,cf$tag,"_sc",sc,subsc,"_cor",cor_val,"_N",target_N,"_alternative", rct_ctrl_tag, ".RData")
          }

          if (file.exists(mbart_file)) {
            mBART_res <- readRDS(mbart_file)
            mBART_res <- mBART_res[mBART_res$alpha == 0.95 & mBART_res$beta == 2, ]
            mBART_res_ATE <- mBART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            mBART_res_sigma <- mBART_res[, c("c","iteration",
                                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                  "pehe_subj","bias_subj")]
            mBART_res_ATE$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
            mBART_res_sigma$Method <- paste0("MAP-BART(", cf$lab, "N=", target_N, ")")
            mBART_res_ATE$SubScenario <- subsc_label
            mBART_res_sigma$SubScenario <- subsc_label
            subsc_res_list_ATE[[length(subsc_res_list_ATE) + 1]] <- mBART_res_ATE
            subsc_res_list_sigma[[length(subsc_res_list_sigma) + 1]] <- mBART_res_sigma
          }
        }, error = function(e) {
          # Silently skip if file doesn't exist
        })
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
    list(name = "LMv1", file_prefix = "LMv1"),
    list(name = "LMv2", file_prefix = "LMv2"),
    list(name = "LMv3", file_prefix = "LMv3"),
    list(name = "HierLM", file_prefix = "HierLM"),
    list(name = "MAP", file_prefix = "MAP"),
    list(name = "PSCL", file_prefix = "PSCL"),
    list(name = "BARTv1", file_prefix = "BARTv1"),
    list(name = "BARTv2", file_prefix = "BARTv2"),
    list(name = "BARTv3", file_prefix = "BARTv3")
  )

  for (subsc in subscenarios) {
    for (cor_val in cor_vals) {
      # Construct null file suffix based on subsc
      if (subsc == "NULL") {
        null_file_suffix <- paste0("_sc", sc, "_cor", cor_val, "_null", rct_ctrl_tag, ".RData")
        subsc_label <- "NULL"
      } else {
        null_file_suffix <- paste0("_sc", sc, subsc, "_cor", cor_val, "_null", rct_ctrl_tag, ".RData")
        subsc_label <- subsc
      }

      for (method_info in methods_list) {
        tryCatch({
          null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
          if (file.exists(null_file)) {
            null_res <- readRDS(null_file)

            # For BART, filter by alpha and beta
            if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3")) {
              null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            }

            # Calculate mean FP (Type I error) - use na.rm = TRUE for PSCL which may have NAs
            if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
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

      # Try to load MAP-BART null results for each control-config x N_target value
      for (cf in .mb_cfgs) for (target_N in mbart_target_Ns) {
        tryCatch({
          if (subsc == "NULL") {
            null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, cf$tag, "_sc", sc, "_cor", cor_val, "_N", target_N, "_null", rct_ctrl_tag, ".RData")
          } else {
            null_file <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, cf$tag, "_sc", sc, subsc, "_cor", cor_val, "_N", target_N, "_null", rct_ctrl_tag, ".RData")
          }
          if (file.exists(null_file)) {
            null_res <- readRDS(null_file)
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                SubScenario = subsc_label,
                c = cor_val,
                Method = paste0("MAP-BART(", cf$lab, "N=", target_N, ")"),
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {
          # Silently skip if null file doesn't exist
        })
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

# Create method levels dynamically.  HierLM still labels by prior value;
# MAP and MAP-BART label by N_target.
lmv2_methods  <- paste0("LMv2(N=", round(lmv2_rwd_w_vals * n_rwd_nominal), ")")
lmv4_methods  <- paste0("HierLM(",      lmv4_prior_vals, ")")
map_methods   <- paste0("MAP(N=",     map_target_Ns,   ")")
mbart_methods <- unlist(lapply(.mb_cfgs, function(cf) paste0("MAP-BART(", cf$lab, "N=", mbart_target_Ns, ")")))
all_method_levels <- c("LMv1", lmv2_methods, "LMv3", "BARTv1", "BARTv2", "BARTv3", "PSCL", map_methods, lmv4_methods, mbart_methods)

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
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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
        filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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

    # Population mean plots
    res_pop_mean_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, map_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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

    summary_table_pop_mean <- res_pop_mean_long %>%
      group_by(SubScenario_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean)

    # Split into treatment and control tables
    summary_table_pop_mean_trt <- summary_table_pop_mean %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_mean_ctrl <- summary_table_pop_mean %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_pop_mean_trt_sep <- add_subsc_separators(summary_table_pop_mean_trt, "SubScenario_Label")
    summary_table_pop_mean_ctrl_sep <- add_subsc_separators(summary_table_pop_mean_ctrl, "SubScenario_Label")
    table_grob_pop_mean_trt <- gridExtra::tableGrob(summary_table_pop_mean_trt_sep, rows = NULL)
    table_grob_pop_mean_ctrl <- gridExtra::tableGrob(summary_table_pop_mean_ctrl_sep, rows = NULL)
    table_grob_pop_mean <- wrap_plots(table_grob_pop_mean_trt, table_grob_pop_mean_ctrl, ncol = 2)

    facet_plot_pop_mean <- ggplot(res_pop_mean_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Mean Performance",
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

    # Combine ATE, population mean, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 1.5, 4, 2))
    } else {
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 2))
    }
  } else {
    # Only ATE plot if no sigma data
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = c(4, 1))
  }

  # Dynamically calculate dimensions based on number of sub-scenarios
  n_subsc <- length(unique(res$SubScenario_Label))
  n_plots <- if (nrow(res_sigma) > 0) 4 else 1  # ATE, subj, pop_mean, sigma or just ATE

  # Calculate total dimensions
  total_height <- n_subsc * 15 * n_plots + 10
  total_width <- 18

  ggsave(paste0(file.path(mainDir),"/mapbart-sim-gaussian-realcov/inserts/p",p_obs,"_sc",sc,rct_ctrl_tag,"_all_subsc_results.jpg"),
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
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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
        filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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

    # Population mean plots
    res_pop_mean_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, map_methods, "BARTv1", "BARTv2", mbart_methods)) %>%
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

    summary_table_pop_mean <- res_pop_mean_long %>%
      group_by(SubSc_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      rename(Correlation = c) %>%
      # Order by subscenario first, then correlation
      arrange(SubSc_Cor_Label, Method)

    summary_table_pop_mean_trt <- summary_table_pop_mean %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_mean_ctrl <- summary_table_pop_mean %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_pop_mean_trt_sep <- add_subsc_separators(summary_table_pop_mean_trt, "SubSc_Cor_Label")
    summary_table_pop_mean_ctrl_sep <- add_subsc_separators(summary_table_pop_mean_ctrl, "SubSc_Cor_Label")
    table_grob_pop_mean_trt <- gridExtra::tableGrob(summary_table_pop_mean_trt_sep, rows = NULL)
    table_grob_pop_mean_ctrl <- gridExtra::tableGrob(summary_table_pop_mean_ctrl_sep, rows = NULL)
    table_grob_pop_mean <- wrap_plots(table_grob_pop_mean_trt, table_grob_pop_mean_ctrl, ncol = 2)

    facet_plot_pop_mean <- ggplot(res_pop_mean_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Mean Performance Across Sub-Scenarios and Correlations",
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

    # Combine ATE, population mean, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 1.5, 4, 2))
    } else {
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = c(4, 2, 4, 2, 4, 2))
    }
  } else {
    # Only ATE plot if no sigma data
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = c(4, 1))
  }

  # Dynamically calculate dimensions based on number of sub-scenarios and correlations
  n_subsc <- length(unique(res$SubSc_Cor_Label))
  n_plots <- if (nrow(res_sigma) > 0) 4 else 1  # ATE, subj, pop_mean, sigma or just ATE

  # Calculate total dimensions
  total_height <- n_subsc * 12 * n_plots + 10
  total_width <- 18

  ggsave(paste0(file.path(mainDir),"/mapbart-sim-gaussian-realcov/inserts/p",p_obs,"_sc",sc,rct_ctrl_tag,"_all_subsc_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}

}, error = function(e) {
  message(sprintf("Skipped sc = %d: %s", sc, conditionMessage(e)))
})
}
