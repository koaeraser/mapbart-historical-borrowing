rm(list=ls())
mainDir <- "/Users/oliviazhang/Desktop/"
library(RBesT)

data_folder <- "data_v2"  # Options: "data" or "data_v2"

sc <- 2
subsc <- "E"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-gaussian", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-gaussian", data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# Calibrated s^2_star comes from ess_local/ess_calibration_map.R, which
# saves res/ess_calibration_map_<version>.RData containing $S_grid and
# $ESS.  For each N_target multiplier m, we pick the smallest s^2 in
# the grid with ESS(s^2) <= m * N_target.  The result is two parallel
# vectors: prior_vals (s^2 values used inside gMAP) and target_Ns (the
# integer target N each was calibrated against -- used in filenames
# and labels in place of the raw s^2 value).
N_target_multipliers <- c(1.00, 0.75, 0.50)

.cal_dir   <- file.path(mainDir, "mapbart-sim-gaussian", "ess_local", "res")
.cal_files <- list.files(.cal_dir,
                           pattern = "^ess_calibration_map_.*\\.RData$",
                           full.names = TRUE)
.cal_files <- .cal_files[!grepl("_checkpoint\\.RData$", .cal_files)]
if (length(.cal_files)) {
  .cal_latest <- .cal_files[order(file.info(.cal_files)$mtime,
                                    decreasing = TRUE)][1]
  .cal        <- readRDS(.cal_latest)
  .pick <- function(m) {
    ok <- .cal$ESS <= m * .cal$N_target
    if (any(ok)) .cal$S_grid[which(ok)[1]] else NA_real_
  }
  prior_vals <- sapply(N_target_multipliers, .pick)
  target_Ns  <- as.integer(round(N_target_multipliers * .cal$N_target))
  ok_idx     <- which(!is.na(prior_vals))
  prior_vals <- prior_vals[ok_idx]
  target_Ns  <- target_Ns [ok_idx]
  cat(sprintf("MAP.R: from %s\n", basename(.cal_latest)))
  for (k in seq_along(prior_vals))
    cat(sprintf("  N=%d -> prior(s^2) = %.4g\n",
                target_Ns[k], prior_vals[k]))
} else {
  prior_vals <- 0.10
  target_Ns  <- 100L
  cat(sprintf("MAP.R: prior_vals = %.4g, target_Ns = %d  (fallback; no calibration file)\n",
              prior_vals, target_Ns))
}

threshold <- 0.95

# Strength of unmeasured confounding
if (sc == 1 | sc == 2){
  cor <- 1
}
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[1]
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
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, subsc, "_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, "_cor", cor, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, subsc, "_cor", cor, "_threshold.RData")
    }
  }

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}

for (.idx in seq_along(prior_vals)) {

  prior    <- prior_vals[.idx]
  target_N <- target_Ns[.idx]

  cat("\n========== Running with N_target =", target_N,
      "  prior(s^2) =", prior, "==========\n")

  set.seed(6)
  seed <- sample(1:10000, niter*4, replace = F)
  ee <- 1

  res <- data.frame(
    dist = "HalfNormal",
    c = cor
  )
  res <- res[rep(1:nrow(res), each = niter), ]
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
  res$pehe_subj <- NA  # Not applicable for MAP
  res$bias_subj <- NA  # Not applicable for MAP
  res$bias.trt.sigma <- NA  # MAP method doesn't estimate these directly
  res$sd.trt.sigma <- NA
  res$w2distance.trt.sigma <- NA
  res$bias.ctrl.sigma <- NA
  res$sd.ctrl.sigma <- NA
  res$w2distance.ctrl.sigma <- NA
  res$bias.trt.median.pop <- NA  # Not applicable for MAP
  res$w2distance.trt.median.pop <- NA
  res$bias.ctrl.median.pop <- NA  # Not applicable for MAP
  res$w2distance.ctrl.median.pop <- NA
  res$iteration <- rep(1:niter, times = nrow(res) / niter)

  # For threshold calibration under H0
  if (hypo == "null") {
    all_decisions <- numeric(niter * length(cor))
    decision_idx <- 1
  }

  cat("Loaded threshold from H0 run:", threshold)

  for (c in cor){
    for (iter in 1:niter){
      
      set.seed(seed[ee])
      ee <- ee + 1
      
      # Construct filename for RCT data (current iteration)
      # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
      if (sc == 1 | sc == 2) {
        if (is.null(subsc)) {
          filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
        } else {
          filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
        }
      } else {
        if (is.null(subsc)) {
          filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
        } else {
          filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
        }
      }
      
      # Construct filename for RWD data (current iteration)
      if (sc == 1 | sc == 2) {
        if (is.null(subsc)) {
          filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
        } else {
          filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
        }
      } else {
        if (is.null(subsc)) {
          filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
        } else {
          filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
        }
      }
      
      # Read RCT data from current iteration
      data_rct_full <- readRDS(filename_rct)
      rct_idx <- data_rct_full$X[,"D"] == 1
      
      # Read RWD data from current iteration
      data_rwd_full <- readRDS(filename_rwd)
      rwd_idx <- data_rwd_full$X[,"D"] == 0
      
      # Combine RCT and RWD from current iteration
      data <- data.frame(
        rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
              data_rwd_full$X[rwd_idx, paste0("X",1:10)])
      )
      data$D <- c(data_rct_full$X[rct_idx, "D"],
                  data_rwd_full$X[rwd_idx, "D"])
      data$Z <- c(data_rct_full$X[rct_idx, "Z"],
                  data_rwd_full$X[rwd_idx, "Z"])
      data$Y <- c(data_rct_full$y[rct_idx],
                  data_rwd_full$y[rwd_idx])
      
      # Also need filename variable for reading true values later
      filename <- filename_rct
      
      print(dim(data))
      
      # Control
      historical_mean <- mean(data[data$D == 0 & data$Z == 0, "Y"])  
      historical_sd   <- sd(data[data$D == 0 & data$Z == 0, "Y"])  
      historical_n    <- nrow(data[data$D == 0 & data$Z == 0, ])
      
      historical_data <- data.frame(
        y = historical_mean, 
        n = historical_n,
        y.se = historical_sd / sqrt(historical_n), 
        study = factor(1:length(historical_mean)) 
      )
      
      map_prior <- gMAP(
        formula = cbind(y, y.se) ~ 1 | study,  
        data = historical_data,
        beta.prior = 10,  
        # tau.dist = "HalfNormal",
        # tau.prior = c(0, 1),
        tau.dist = "InvGamma",
        tau.prior = c(3/2, 3*prior/2),
        family = gaussian
      )
      
      # samples <- unlist(map_prior$fit@sim$samples[[1]]["theta_pred"])
      map_mixture <- automixfit(map_prior)
      print(ess(map_mixture, sigma = 1.5))
      
      new_mean <- mean(data[data$D == 1 & data$Z == 0, "Y"])  
      new_sd   <- sd(data[data$D == 1 & data$Z == 0, "Y"])  
      new_n    <- nrow(data[data$D == 1 & data$Z == 0, ])
      
      posterior_ctrl <- postmix(map_mixture, 
                                m = new_mean, 
                                se = new_sd / sqrt(new_n))
      
      # Treatment
      treatment_mean <- mean(data[data$D == 1 & data$Z == 1, "Y"])  
      treatment_sd   <- sd(data[data$D == 1 & data$Z == 1, "Y"])  
      treatment_n    <- nrow(data[data$D == 1 & data$Z == 1, ])
      
      # Define a weakly informative normal prior for treatment
      weak_prior <- mixnorm(c(1,              #w
                              0,              #mean
                              10              #se
      ), param = 'ms')
      posterior_trt <- postmix(weak_prior, 
                               m = treatment_mean, 
                               se = treatment_sd / sqrt(treatment_n))
      
      #-------------------------------------
      #---------- Calculate ATE ------------
      #-------------------------------------
      samples_trt <- rmix(posterior_trt, 1500)
      samples_ctrl <- rmix(posterior_ctrl, 1500)
      
      post_samples <- samples_trt - samples_ctrl
      eff <- readRDS(filename)$treat_eff
      eff_star <- readRDS(filename)$treat_eff_star
      
      # Posterior estimate
      delta_hat <- mean(post_samples, na.rm = TRUE)
      
      # Within replicate metrics
      bias <- delta_hat - eff
      w2distance <- sqrt(mean((post_samples - eff)^2, na.rm = TRUE))
      
      # Across replicate metrics
      rmse <- (delta_hat - eff)^2
      
      # Other
      SD <- sqrt(var(post_samples, na.rm = TRUE))
      w1distance <- mean(abs(post_samples - eff), na.rm = TRUE)
      coverage <- quantile(post_samples, probs = c(0.025), na.rm = TRUE) <= eff &
        quantile(post_samples, probs = c(0.975), na.rm = TRUE) >= eff
      ci <- diff(quantile(post_samples, c(0.025, 0.975), na.rm = TRUE))
      
      # True Positive/Power
      decision <- mean(post_samples > eff_star, na.rm = TRUE)
      
      tp_calibrated <- as.numeric(decision > threshold)
      tp <- as.numeric(decision > 0.95)
      
      if (hypo == "null") {
        all_decisions[decision_idx] <- decision
        decision_idx <- decision_idx + 1
      }
      
      # False Positive
      if (hypo == "null") {
        fp <- tp
      } else {
        fp <- NA
      }
      
      res[which(res$c==c & res$iteration==iter), "bias"] <- bias
      res[which(res$c==c & res$iteration==iter), "sd"] <- SD
      res[which(res$c==c & res$iteration==iter), "rmse"] <- rmse
      res[which(res$c==c & res$iteration==iter), "w1distance"] <- w1distance
      res[which(res$c==c & res$iteration==iter), "w2distance"] <- w2distance
      res[which(res$c==c & res$iteration==iter), "ci"] <- ci
      res[which(res$c==c & res$iteration==iter), "coverage"] <- coverage
      res[which(res$c==c & res$iteration==iter), "tp_calibrated"] <- tp_calibrated
      res[which(res$c==c & res$iteration==iter), "tp"] <- tp
      res[which(res$c==c & res$iteration==iter), "fp"] <- fp
      
      # Subject-wise calculations - not applicable for MAP (operates at group level)
      res[which(res$c==c & res$iteration==iter), "pehe_subj"] <- NA
      res[which(res$c==c & res$iteration==iter), "bias_subj"] <- NA
      
      # Variance estimation for treatment and control groups
      # MAP method doesn't model residual variance, so these remain NA
      res[which(res$c==c & res$iteration==iter), "bias.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==iter), "sd.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==iter), "w2distance.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==iter), "bias.ctrl.sigma"] <- NA
      res[which(res$c==c & res$iteration==iter), "sd.ctrl.sigma"] <- NA
      res[which(res$c==c & res$iteration==iter), "w2distance.ctrl.sigma"] <- NA
      
      #--------------------------------------------------
      #------ Population mean estimation ----------------
      #--------------------------------------------------
      true_mean_trt_pop <- readRDS(filename)$true_mean_trt_pop
      true_mean_ctrl_pop <- readRDS(filename)$true_mean_ctrl_pop
      
      # Calculate population means
      mean_trt_pop_hat <- mean(samples_trt)
      mean_ctrl_pop_hat <- mean(samples_ctrl)
      
      # Calculate metrics
      bias_trt_median_pop <- mean_trt_pop_hat - true_mean_trt_pop
      w2distance_trt_median_pop <- sqrt(mean((samples_trt - true_mean_trt_pop)^2))
      bias_ctrl_median_pop <- mean_ctrl_pop_hat - true_mean_ctrl_pop
      w2distance_ctrl_median_pop <- sqrt(mean((samples_ctrl - true_mean_ctrl_pop)^2))
      
      res[which(res$c==c & res$iteration==iter), "bias.trt.median.pop"] <- bias_trt_median_pop
      res[which(res$c==c & res$iteration==iter), "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
      res[which(res$c==c & res$iteration==iter), "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
      res[which(res$c==c & res$iteration==iter), "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop
      
      print(paste0("Done for iteration ", iter, " cor = ", c))
    }
    print(paste0("Done for cor = ", c))
  }

  # Threshold calibration under null hypothesis
  if (hypo == "null") {
    cat("\n=== Threshold Calibration (H0) ===\n")
    cat("Original threshold: ", threshold, "\n")

    # Grid search for exact type I error = 0.05 (more precise)
    candidate_thresholds <- seq(0.1, 0.999, by = 0.001)
    type1_errors <- sapply(candidate_thresholds, function(thresh) {
      mean(all_decisions > thresh, na.rm = TRUE)
    })

    best_idx <- which.min(abs(type1_errors - 0.05))
    calibrated_threshold_grid <- candidate_thresholds[best_idx]
    actual_type1_error <- type1_errors[best_idx]

    cat("Selected threshold: ", calibrated_threshold_grid, "\n")
    cat("Actual Type I error: ", actual_type1_error, "\n")

    # Save calibrated threshold for use in alternative hypothesis run
    # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
    if (sc == 1 | sc == 2) {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, "_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, subsc, "_threshold.RData")
      }
    } else {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, "_cor", c, "_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/MAP_p", p_obs, "_sc", sc, subsc, "_cor", c, "_threshold.RData")
      }
    }
    saveRDS(calibrated_threshold_grid, file = threshold_file)
  }

  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  # Save results with prior value in filename
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,"_N",target_N,"_",hypo,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,subsc,"_N",target_N,"_",hypo,".RData"))
    }
  } else {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,"_cor",c,"_N",target_N,"_",hypo,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,subsc,"_cor",c,"_N",target_N,"_",hypo,".RData"))
    }
  }

} # End of (prior, target_N) loop
