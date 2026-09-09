rm(list=ls())
# devtools::install_github("olssol/psrwe")
library(miceadds)
library(dplyr)
mainDir <- "/Users/oliviazhang/Desktop/"
data_folder <- "data_v2"  # Options: "data" or "data_v2"

old_wd <- getwd()
setwd(paste0(mainDir, "psrwe"))
source.all("R/")
setwd(old_wd)

sc <- 3
subsc <- "E"  # NULL, "A", "B", "C", "E", or "F"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-gaussian", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-gaussian", data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

threshold <- 0.95

# Strength of unmeasured confounding
if (sc == 1 | sc == 2){
  cor <- 1
}
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[3]
}

# Infer niter from the number of simulated data files matching this scenario
data_prefix <- if (sc == 1 | sc == 2) {
  paste0("data_p", p_obs, "_sc", sc, if (!is.null(subsc)) subsc else "", "_", hypo, "_")
} else {
  paste0("data_p", p_obs, "_sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor, "_", hypo, "_")
}
niter <- length(Sys.glob(file.path(mainDir, "mapbart-sim-gaussian", data_folder, paste0(data_prefix, "*.RData"))))
if (niter == 0) stop("No simulated data files found for this scenario")

# Auto-load calibrated threshold when running under alternative
if (hypo == "alternative") {
  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/PSCL_p", p_obs, "_sc", sc, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/PSCL_p", p_obs, "_sc", sc, subsc, "_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/PSCL_p", p_obs, "_sc", sc, "_",
                               cor, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/PSCL_p", p_obs, "_sc", sc, subsc, "_",
                               cor, "_threshold.RData")
    }
  }

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}

set.seed(6)
seed <- sample(1:10000, niter*4, replace = F)
ee <- 1

res <- data.frame(c = rep(cor, each = niter))
res$bias <- NA
res$sd <- NA
res$rmse <- NA
res$w1distance <- NA
res$w2distance <- NA
res$ci <- NA
res$coverage <- NA
res$tp_calibrated <- NA
res$tp <- NA
res$fp <- NA
res$pehe_subj <- NA  # Not applicable for PSCL
res$bias_subj <- NA  # Not applicable for PSCL
res$bias.trt.sigma <- NA  # Not applicable for PSCL
res$sd.trt.sigma <- NA
res$w2distance.trt.sigma <- NA
res$bias.ctrl.sigma <- NA
res$sd.ctrl.sigma <- NA
res$w2distance.ctrl.sigma <- NA
res$bias.trt.median.pop <- NA  # Not applicable for PSCL
res$w2distance.trt.median.pop <- NA
res$bias.ctrl.median.pop <- NA  # Not applicable for PSCL
res$w2distance.ctrl.median.pop <- NA
res$iteration <- rep(1:niter, times = nrow(res) / niter)

# For threshold calibration under H0
if (hypo == "null") {
  all_decisions <- numeric(niter * length(cor))
  decision_idx <- 1
}

cat("Loaded threshold from H0 run:", threshold)

for (c_idx in seq_along(cor)){
  c <- cor[c_idx]
for (iter in 1:niter){

  set.seed(seed[ee])
  ee <- ee + 1

  # Calculate row index directly
  row_idx <- (c_idx - 1) * niter + iter

  # Use tryCatch to compute results, return list on success or NULL on error
  iter_result <- tryCatch({

    # Construct filename for RCT data (current iteration)
    # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
    if (sc == 1 | sc == 2) {
      if (!is.null(subsc)) {
        filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
      } else {
        filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
      }
    } else {
      if (!is.null(subsc)) {
        filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
      } else {
        filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
      }
    }

    # Construct filename for RWD data (current iteration)
    if (sc == 1 | sc == 2) {
      if (!is.null(subsc)) {
        filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
      } else {
        filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
      }
    } else {
      if (!is.null(subsc)) {
        filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
      } else {
        filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
      }
    }

    # Read RCT data from current iteration
    data_rct_full <- readRDS(filename_rct)
    rct_idx <- data_rct_full$X[,"D"] == 1

    # Read RWD data from current iteration
    data_rwd_full <- readRDS(filename_rwd)
    rwd_idx <- data_rwd_full$X[,"D"] == 0

    # Combine RCT and RWD from current iteration
    data <- cbind(
      rbind(data_rct_full$X[rct_idx, ],
            data_rwd_full$X[rwd_idx, ]),
      data.frame("y" = c(data_rct_full$y[rct_idx],
                         data_rwd_full$y[rwd_idx]))
    )

    # Also need filename variable for reading true values later
    filename <- filename_rct

    ps <- psrwe_est(data, v_covs = paste("X", 1:10, sep = ""),
                    v_grp = "D", cur_grp_level = "1",
                    v_arm = "Z", ctl_arm_level = "0")

    borrow <- psrwe_borrow(ps, total_borrow = 100,
                           method = "distance")
    fit <- psrwe_compl(borrow, v_outcome = "y",
                       outcome_type = "continuous")
    fit_ci <- psrwe_outana(fit)

    #-------------------------------------
    #---------- Calculate ATE ------------
    #-------------------------------------
    eff <- readRDS(filename)$treat_eff
    eff_star <- readRDS(filename)$treat_eff_star

    # Posterior estimate (point estimate for PSCL)
    delta_hat <- fit$Effect$Overall_Estimate$Mean
    print(delta_hat)
    
    # Within replicate metrics
    bias <- delta_hat - eff
    w2distance <- NA

    # Across replicate metrics
    rmse <- (delta_hat - eff)^2

    # Other
    SD <- fit$Effect$Overall_Estimate$StdErr
    w1distance <- NA
    ci_lower <- fit_ci$CI$Effect$Overall_Estimate$Lower
    ci_upper <- fit_ci$CI$Effect$Overall_Estimate$Upper
    ci <- ci_upper - ci_lower
    coverage <- ci_lower <= eff & ci_upper >= eff

    if (hypo == "null") {
      fp <- as.numeric(ci_lower > eff_star)
      tp <- NA
      tp_calibrated <- NA
    } else {
      tp <- as.numeric(ci_lower > eff_star)
      fp <- NA
      tp_calibrated <- NA
    }

    # Return results as a list
    list(bias = bias, sd = SD, rmse = rmse, w1distance = w1distance,
         w2distance = w2distance, ci = ci, coverage = coverage,
         tp_calibrated = tp_calibrated, tp = tp, fp = fp)

  }, error = function(e) {
    cat(paste0("\nError in iteration ", iter, " for cor = ", c, ":\n"))
    cat(paste0("  ", e$message, "\n"))
    cat("  Skipping this iteration and setting results to NA\n\n")
    NULL  # Return NULL on error
  })

  # Assign results outside tryCatch (if successful)
  if (!is.null(iter_result)) {
    res[row_idx, "bias"] <- iter_result$bias
    res[row_idx, "sd"] <- iter_result$sd
    res[row_idx, "rmse"] <- iter_result$rmse
    res[row_idx, "w1distance"] <- iter_result$w1distance
    res[row_idx, "w2distance"] <- iter_result$w2distance
    res[row_idx, "ci"] <- iter_result$ci
    res[row_idx, "coverage"] <- iter_result$coverage
    res[row_idx, "tp_calibrated"] <- iter_result$tp_calibrated
    res[row_idx, "tp"] <- iter_result$tp
    res[row_idx, "fp"] <- iter_result$fp
  }
  # If iter_result is NULL (error), res stays NA (default)

  print(paste0("Done for iteration ", iter, " cor = ", c))
}
  print(paste0("Done for cor = ", c))
}

# Note: PSCL method doesn't support threshold calibration as it doesn't provide posterior samples
# No threshold calibration section needed for PSCL

# Note: RMSE is already set to NA for each iteration above since PSCL
# doesn't provide posterior samples needed to calculate it properly

# For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
if (sc == 1 | sc == 2) {
  if (is.null(subsc)) {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/PSCL_p",p_obs,"_sc",sc,"_",hypo,".RData"))
  } else {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/PSCL_p",p_obs,"_sc",sc,subsc,"_",hypo,".RData"))
  }
} else {
  if (is.null(subsc)) {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/PSCL_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,".RData"))
  } else {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/PSCL_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,".RData"))
  }
}
