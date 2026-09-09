rm(list = ls())
library(MASS)
library(Rlab)

mainDir <- "/Users/oliviazhang/Desktop/"

sc <- 3
subsc <- "E"
hypo <- "alternative"

if (hypo == "null"){
  niter <- 500
} else{
  niter <- 500
}
n <- 300

r <- 2
n1 <- (r/(r+1)) * n
n3 <- n

n_rct_pool <- 2 * (n1 + n3)
n_rwd_pool <- 2 * (n1 + n3)

p_obs <- 10
conf_cols <- c(5, 6)
extreme_cols <- c(1)

if (sc == 1 | sc == 2){
  cor <- 1
}
if (sc == 3){
  cor <- 0.7
}

set.seed(6)
seed <- sample(1:10000, niter*length(cor), replace = F)
ee <- 1

seed_rwd <- 6

# RWD-generation mode: TRUE = one frozen pool (files tagged _fz); FALSE = fresh RWD each iteration
rwd_frozen <- TRUE
set.seed(seed_rwd)
seed_rwd_iter <- sample(1:10000, niter*length(cor), replace = F)   # per-iteration RWD seeds

sd_x <- 1
sd_U <- 1
m_x <- c(1, -1, 0.5, -1, 2, -2, 2, -3, 1, 0)

if (sc == 2) {
  beta_D   <- -0.55
  sel_temp <- 1.0
}
if (sc == 3) {
  beta_D   <- -0.55
  sel_temp <- 1.0
}
if (hypo == "null") {
  eff <- 0
} else {
  eff <- 1.0
}
eff_star <- 0.5

if (subsc == "E"){
  if (sc == 1) {
    b0_rct <- 0.4
    b0_rwd <- 0.4

    sd_Y_rct <- 1.5
    sd_Y_rwd <- 1.5

  } else if (sc == 2) {
    b0_rct <- -0.1
    b0_rwd <- -0.1

    sd_Y_rct <- 1.5
    sd_Y_rwd <- 1.5

  } else if (sc == 3) {
    b0_rct <- -0.1
    b0_rwd <- -0.1

    sd_Y_rct <- 1.5
    sd_Y_rwd <- 1.5

  }

  beta_rct <- c(-0.50, -0.75, -0.50, -0.50,
                -1.35, -0.80, -0.01, -0.01, -0.01, -0.01)
  beta_rwd <- beta_rct

  gamma <- 0
  eta <- 0

}

#===========================================================================
#============================ Generate RWD_pool ============================
#===========================================================================
# Factored so the RWD pool can be built once (frozen) or refreshed per iteration.
# Assigns the pool objects (X_rwd_pool/W_rwd_pool/y_rwd_pool/lp_rwd_pool) to the global scope.
gen_rwd_pool <- function(seed_x, seed_y) {
  #========================== Generate X ==========================
  set.seed(seed_x)
  Xp <- matrix(rep(NA, n_rwd_pool * p_obs), ncol = p_obs)
  if (sc == 3) {
    Up <- matrix(rep(NA, n_rwd_pool * length(conf_cols)), ncol = length(conf_cols))
  }
  for (i in 1:n_rwd_pool) {
    Xp[i, 1:p_obs] <- rnorm(p_obs, mean = m_x[1:p_obs], sd = sd_x)
    for (col in extreme_cols) {
      if (runif(1) < 0.15) {
        Xp[i, col] <- rnorm(1, mean = sample(c(-3, 3), 1), sd = 0.5)
      }
    }

    #========================== Generate U ==========================
    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        Up[i, j] <- rnorm(1, mean = m_x[conf_cols[j]], sd = sd_x)
      }
    }
  }

  #========================== Generate Y ==========================
  set.seed(seed_y)
  if (sc == 3) {
    Wp <- Up
  } else {
    Wp <- Xp[, conf_cols]
  }
  yp <- rep(NA, n_rwd_pool)
  lpp <- rep(NA, n_rwd_pool)
  for (i in 1:n_rwd_pool) {
    x_vec <- Xp[i, ]
    x_vec[conf_cols] <- Wp[i, ]

    lp <- b0_rwd +
      beta_rwd[1] * x_vec[1]^2 +
      beta_rwd[2] * exp(-(x_vec[2] + 1)^2 / 2) +
      beta_rwd[3] * abs(x_vec[3] - 0.5) +
      beta_rwd[4] * abs(x_vec[4] + 1) +
      beta_rwd[5] * tanh(x_vec[5] - 2) +
      beta_rwd[6] * tanh(x_vec[6] + 2) +
      2.0 * ( (x_vec[5] - 2)^2 / (1 + ((x_vec[5] - 2)/1.5)^2) +
              (x_vec[6] + 2)^2 / (1 + ((x_vec[6] + 2)/1.5)^2) ) +
      beta_rwd[7] * x_vec[7] +
      beta_rwd[8] * x_vec[8] +
      beta_rwd[9] * x_vec[9] +
      beta_rwd[10] * x_vec[10]

    lpp[i] <- lp
    yp[i] <- lp + rnorm(n = 1, sd = sd_Y_rwd)
  }

  # Build the observed X pool for each correlation level (sc == 3): the
  # conf_cols of the observed X are proxies correlated with the latent U.
  X_rwd_pool_by_c <- vector("list", length(cor))
  for (ci in seq_along(cor)) {
    Xc <- Xp
    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        U_std <- (Up[,j] - mean(Up[,j])) / sd(Up[,j])
        set.seed(seed_x + j + 100)
        Xc[, conf_cols[j]] <- cor[ci] * U_std + sqrt(1 - cor[ci]^2) * rnorm(n_rwd_pool, mean = 0, sd = sd_U)
      }
    }
    X_rwd_pool_by_c[[ci]] <- Xc
  }

  X_rwd_pool      <<- Xp
  W_rwd_pool      <<- Wp
  y_rwd_pool      <<- yp
  lp_rwd_pool     <<- lpp
  X_rwd_pool_by_c <<- X_rwd_pool_by_c
  if (sc == 3) U_rwd_pool <<- Up
}

# Frozen mode: build the pool once, up front (per-iteration mode rebuilds it inside the loop).
if (rwd_frozen) gen_rwd_pool(seed_rwd, seed_rwd + 1000)

for (ci in seq_along(cor)) {
  c <- cor[ci]
  if (rwd_frozen) X_rwd_pool <- X_rwd_pool_by_c[[ci]]

  #===========================================================================
  #============================ Generate RCT_pool ============================
  #===========================================================================
  for (iter in 1:niter){

    # Per-iteration mode: fresh RWD pool each iteration. Kept BEFORE set.seed(seed[ee])
    # so the RCT stream (and thus the RCT data) is identical to frozen mode.
    if (!rwd_frozen) { gen_rwd_pool(seed_rwd_iter[ee], seed_rwd_iter[ee] + 1000); X_rwd_pool <- X_rwd_pool_by_c[[ci]] }

    set.seed(seed[ee])
    ee <- ee + 1

    #========================== Generate X ==========================
    X_rct_pool <- matrix(rep(NA, n_rct_pool * p_obs), ncol = p_obs)
    if (sc == 3) {
      U_rct_pool <- matrix(rep(NA, n_rct_pool * length(conf_cols)), ncol = length(conf_cols))
    }

    for (i in 1:n_rct_pool) {
      X_rct_pool[i, 1:p_obs] <- rnorm(p_obs, mean = m_x[1:p_obs], sd = sd_x)

      for (col in extreme_cols) {
        if (runif(1) < 0.15) {
          X_rct_pool[i, col] <- rnorm(1, mean = sample(c(-3, 3), 1), sd = 0.5)
        }
      }

      if (sc == 3) {
        for (j in 1:length(conf_cols)) {
          U_rct_pool[i, j] <- rnorm(1, mean = m_x[conf_cols[j]], sd = sd_x)
        }
      }

    }

    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        U_std <- (U_rct_pool[,j] - mean(U_rct_pool[,j])) / sd(U_rct_pool[,j])
        X_rct_pool[, conf_cols[j]] <- c * U_std + sqrt(1 - c^2) * rnorm(n_rct_pool, mean = 0, sd = sd_U)
      }
    }

    # Combine RCT and RWD
    X <- rbind(X_rct_pool, X_rwd_pool)
    if (sc == 3) {
      U <- rbind(U_rct_pool, U_rwd_pool)
    }

    #========================== Generate D ==========================
    n_total_pool <- n_rct_pool + n_rwd_pool
    if (sc == 1) {
      D <- rbern(n_total_pool, prob = n1/(n1+n3))
    }
    if (sc == 2) {
      if (length(beta_D) > 1) {
        lp <- X[, conf_cols] %*% beta_D
      } else {
        lp <- beta_D * rowSums(X[, conf_cols])
      }
      thres <- quantile(lp, n1/(n1+n3))
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }
    if (sc == 3) {
      if (length(beta_D) > 1) {
        lp <- U[, 1:length(conf_cols)] %*% beta_D
      } else {
        lp <- beta_D * rowSums(U[, 1:length(conf_cols)])
      }
      thres <- quantile(lp, n1/(n1+n3))
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }
    
    #===========================================================================
    #============================== Generate RCT ===============================
    #===========================================================================

    idx_rct_selected <- which(D == 1 & (1:n_total_pool) <= n_rct_pool)
    
    if (length(idx_rct_selected) >= n1) {
      idx_rct <- idx_rct_selected[1:n1]
    } else {
      idx_rct <- idx_rct_selected
    }
    
    n_rct <- length(idx_rct)
    X_rct_sel <- X[idx_rct, , drop = FALSE]
    if (sc == 3) {
      W_rct_sel <- matrix(U[idx_rct, ], ncol = length(conf_cols))
    } else {
      W_rct_sel <- X_rct_sel[, conf_cols]
    }
    
    #========================== Generate Z ==========================
    Z_rct <- rep(1, n_rct)
    
    #========================== Generate Y ==========================
    y_rct <- rep(NA, n_rct)
    lp_rct <- rep(NA, n_rct)
    
    true_mean_trt <- rep(NA, n_rct)
    true_mean_ctrl <- rep(NA, n_rct)
    eff_i_rct <- rep(NA, n_rct)
    
    for (i in 1:n_rct) {
      x_vec <- X_rct_sel[i, ]
      x_vec[conf_cols] <- W_rct_sel[i, ]
      modifier <- gamma * sum(x_vec[conf_cols])
      
      lp <- b0_rct +
        beta_rct[1] * x_vec[1]^2 +
        beta_rct[2] * exp(-(x_vec[2] + 1)^2 / 2) +
        beta_rct[3] * abs(x_vec[3] - 0.5) +
        beta_rct[4] * abs(x_vec[4] + 1) +
        beta_rct[5] * tanh(x_vec[5] - 2) +
        beta_rct[6] * tanh(x_vec[6] + 2) +
        2.0 * ( (x_vec[5] - 2)^2 / (1 + ((x_vec[5] - 2)/1.5)^2) +
                (x_vec[6] + 2)^2 / (1 + ((x_vec[6] + 2)/1.5)^2) ) +
        beta_rct[7] * x_vec[7] +
        beta_rct[8] * x_vec[8] +
        beta_rct[9] * x_vec[9] +
        beta_rct[10] * x_vec[10]
      
      eff_i <- eff
      
      if(Z_rct[i] == 1){
        lp_rct[i] <- lp + modifier + eff_i
        y_rct[i] <- lp + modifier + eff_i + rnorm(n = 1, sd = sd_Y_rct)
      } else {
        lp_rct[i] <- lp + modifier
        y_rct[i] <- lp + modifier + rnorm(n = 1, sd = sd_Y_rct)
      }
      
      eff_i_rct[i] <- eff_i
      true_mean_trt[i] <- lp + modifier + eff_i
      true_mean_ctrl[i] <- lp + modifier
    }
    
    mean_trt_pop <- mean(true_mean_trt)
    mean_ctrl_pop <- mean(true_mean_ctrl)
    
    #===========================================================================
    #============================== Generate RWD ===============================
    #===========================================================================
    idx_rwd_selected <- which(D == 0 & (1:n_total_pool) > n_rct_pool)
    pool_indices_available <- idx_rwd_selected - n_rct_pool
    
    if (length(pool_indices_available) >= n3) {
      pool_indices <- pool_indices_available[1:n3]
    } else {
      pool_indices <- pool_indices_available
    }
    
    n_rwd <- length(pool_indices)
    X_rwd <- X_rwd_pool[pool_indices, , drop = FALSE]
    if (sc == 3) {
      U_rwd <- U_rwd_pool[pool_indices, , drop = FALSE]
    }

    y_rwd <- y_rwd_pool[pool_indices]
    lp_rwd <- lp_rwd_pool[pool_indices]
    Z_rwd <- rep(0, n_rwd)
    
    n_rct_trt <- sum(Z_rct == 1)
    n_rct_ctrl <- sum(Z_rct == 0)
    cat("n_rct_trt =", n_rct_trt, ", n_rct_ctrl =", n_rct_ctrl, ", n_rwd =", n_rwd,"\n")
    
    X_rct_df <- data.frame(X_rct_sel[, 1:p_obs])
    colnames(X_rct_df) <- paste0("X", 1:p_obs)
    X_rwd_df <- data.frame(X_rwd[, 1:p_obs])
    colnames(X_rwd_df) <- paste0("X", 1:p_obs)
    X_rct_df$D <- 1
    X_rwd_df$D <- 0
    X_rct_df$Z <- Z_rct
    X_rwd_df$Z <- Z_rwd

    if (sc == 3) {
      U_rct_df <- data.frame(W_rct_sel)
      colnames(U_rct_df) <- paste0("U", 1:length(conf_cols))
      U_rwd_df <- data.frame(U_rwd)
      colnames(U_rwd_df) <- paste0("U", 1:length(conf_cols))
      U_combined <- rbind(U_rct_df, U_rwd_df)
    } else {
      U_combined <- NULL
    }

    cat("eff =", eff,
        "mean(eff_i_rct) =", mean(eff_i_rct),
        "eff_true =", mean_trt_pop - mean_ctrl_pop,"\n")
    
    # Output
    out <- list(X = rbind(X_rct_df, X_rwd_df),
                U = U_combined,
                y = c(y_rct, y_rwd),
                treat_eff = eff,
                treat_eff_true = mean_trt_pop - mean_ctrl_pop,
                treat_eff_star = eff_star,
                gamma_eff = gamma,  # modifier coefficient
                eta_eff = if (exists("eta")) eta else 0,  # HTE coefficient for treatment effect
                sigma_rct = sd_Y_rct,
                sigma_rwd = sd_Y_rwd,
                true_mean_trt = true_mean_trt,
                true_mean_ctrl = true_mean_ctrl,
                true_mean_trt_pop = mean_trt_pop,
                true_mean_ctrl_pop = mean_ctrl_pop,
                eff_i = eff_i_rct,
                lp = c(lp_rct, lp_rwd),
                n_rct = n_rct,
                n_rwd = n_rwd
    )
    
    # Save with scenario identifiers
    if (sc == 1 | sc == 2) {
      if (is.null(subsc)) {
        saveRDS(out, file = paste0(mainDir, "/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc,"_", hypo, if (rwd_frozen) "_fz" else "", "_", iter, ".RData"))
      } else {
        saveRDS(out, file = paste0(mainDir, "/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc,subsc,"_", hypo, if (rwd_frozen) "_fz" else "", "_", iter, ".RData"))
      }
    } else {
      if (is.null(subsc)) {
        saveRDS(out, file = paste0(mainDir, "/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc,"_cor", c, "_", hypo, if (rwd_frozen) "_fz" else "", "_", iter, ".RData"))
      } else {
        saveRDS(out, file = paste0(mainDir, "/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc,subsc,"_cor", c, "_", hypo, if (rwd_frozen) "_fz" else "", "_", iter, ".RData"))
      }
    }
  }
}

## Testing code - update sc and subsc as needed
sc_test <- sc
subsc_test <- subsc

cat("\n========================================\n")
cat("Testing R² for sc =", sc_test, ", subsc =", ifelse(is.null(subsc_test), "NULL", subsc_test), "\n")
cat("========================================\n")

for (c_test in cor) {

  R2_rct_trt_all <- numeric(niter)
  R2_rct_ctrl_all <- numeric(niter)
  R2_rwd_ctrl_all <- numeric(niter)

  R2_rct_trt_X_all <- numeric(niter)
  R2_rct_ctrl_X_all <- numeric(niter)
  R2_rwd_ctrl_X_all <- numeric(niter)

  # Overlap measures (SMD-based) for X
  smd_trt_ctrl_all <- matrix(NA, nrow = niter, ncol = p_obs)  # RCT trt vs RCT ctrl
  smd_trt_rwd_all <- matrix(NA, nrow = niter, ncol = p_obs)   # RCT trt vs RWD ctrl
  smd_ctrl_rwd_all <- matrix(NA, nrow = niter, ncol = p_obs)  # RCT ctrl vs RWD ctrl

  # Overlap measures (SMD-based) for U (sc == 3 only)
  n_U <- length(conf_cols)
  smd_U_trt_ctrl_all <- matrix(NA, nrow = niter, ncol = n_U)  # RCT trt vs RCT ctrl
  smd_U_trt_rwd_all <- matrix(NA, nrow = niter, ncol = n_U)   # RCT trt vs RWD ctrl
  smd_U_ctrl_rwd_all <- matrix(NA, nrow = niter, ncol = n_U)  # RCT ctrl vs RWD ctrl


  for (i in 1:niter){
    if (sc_test == 1 | sc_test == 2) {
      if (is.null(subsc_test)) {
        filepath <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      } else {
        filepath <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc_test,subsc_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      }
    } else {
      if (is.null(subsc_test)) {
        filepath <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc_test,"_cor",c_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      } else {
        filepath <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/data_v3/data_p",p_obs,"_sc",sc_test,subsc_test,"_cor",c_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      }
    }

    # Read data
    data_sample <- readRDS(filepath)
    X_all <- data_sample$X
    D <- X_all[,"D"]
    Z <- X_all[,"Z"]

    # Use saved lp
    lp_all <- data_sample$lp

    # Calculate R² = Var(lp) / (Var(lp) + sigma^2)
    # For RCT Treatment (D=1, Z=1)
    idx_rct_trt <- D == 1 & Z == 1
    lp_rct_trt <- lp_all[idx_rct_trt]
    R2_rct_trt_all[i] <- var(lp_rct_trt) / (var(lp_rct_trt) + data_sample$sigma_rct^2)

    # For RCT Control (D=1, Z=0)
    idx_rct_ctrl <- D == 1 & Z == 0
    lp_rct_ctrl <- lp_all[idx_rct_ctrl]
    R2_rct_ctrl_all[i] <- var(lp_rct_ctrl) / (var(lp_rct_ctrl) + data_sample$sigma_rct^2)

    # For RWD Control (D=0, Z=0)
    idx_rwd_ctrl <- D == 0 & Z == 0
    lp_rwd_ctrl <- lp_all[idx_rwd_ctrl]
    R2_rwd_ctrl_all[i] <- var(lp_rwd_ctrl) / (var(lp_rwd_ctrl) + data_sample$sigma_rwd^2)

    # Calculate R² using X (observed covariates), sc == 3 only
    if (!is.null(data_sample$U)) {
      X_mat <- as.matrix(X_all[, paste0("X", 1:p_obs)])

      # Calculate lp using X for RCT subjects
      lp_X_rct <- b0_rct +
        beta_rct[1] * X_mat[,1]^2 +
        beta_rct[2] * exp(-(X_mat[,2] + 1)^2 / 2) +
        beta_rct[3] * abs(X_mat[,3] - 0.5) +
        beta_rct[4] * abs(X_mat[,4] + 1) +
        beta_rct[5] * tanh(X_mat[,5] - 2) +
        beta_rct[6] * tanh(X_mat[,6] + 2) +
        2.0 * ( (X_mat[,5] - 2)^2 / (1 + ((X_mat[,5] - 2)/1.5)^2) +
                  (X_mat[,6] + 2)^2 / (1 + ((X_mat[,6] + 2)/1.5)^2) ) +
        beta_rct[7] * X_mat[,7] +
        beta_rct[8] * X_mat[,8] +
        beta_rct[9] * X_mat[,9] +
        beta_rct[10] * X_mat[,10]

      # Calculate lp using X for RWD subjects (use b0_rwd and beta_rwd)
      lp_X_rwd <- b0_rwd +
        beta_rwd[1] * X_mat[,1]^2 +
        beta_rwd[2] * exp(-(X_mat[,2] + 1)^2 / 2) +
        beta_rwd[3] * abs(X_mat[,3] - 0.5) +
        beta_rwd[4] * abs(X_mat[,4] + 1) +
        beta_rwd[5] * tanh(X_mat[,5] - 2) +
        beta_rwd[6] * tanh(X_mat[,6] + 2) +
        2.0 * ( (X_mat[,5] - 2)^2 / (1 + ((X_mat[,5] - 2)/1.5)^2) +
                  (X_mat[,6] + 2)^2 / (1 + ((X_mat[,6] + 2)/1.5)^2) ) +
        beta_rwd[7] * X_mat[,7] +
        beta_rwd[8] * X_mat[,8] +
        beta_rwd[9] * X_mat[,9] +
        beta_rwd[10] * X_mat[,10]

      Y_all <- data_sample$y

      eff_i_X <- eff

      lp_X_rct_trt <- lp_X_rct[idx_rct_trt] + gamma * rowSums(X_mat[idx_rct_trt, conf_cols]) + eff_i_X
      Y_rct_trt <- Y_all[idx_rct_trt]
      R2_rct_trt_X_all[i] <- cor(Y_rct_trt, lp_X_rct_trt)^2

      lp_X_rct_ctrl <- lp_X_rct[idx_rct_ctrl] + gamma * rowSums(X_mat[idx_rct_ctrl, conf_cols])
      Y_rct_ctrl <- Y_all[idx_rct_ctrl]
      R2_rct_ctrl_X_all[i] <- cor(Y_rct_ctrl, lp_X_rct_ctrl)^2

      lp_X_rwd_ctrl <- lp_X_rwd[idx_rwd_ctrl]
      Y_rwd_ctrl <- Y_all[idx_rwd_ctrl]
      R2_rwd_ctrl_X_all[i] <- cor(Y_rwd_ctrl, lp_X_rwd_ctrl)^2
    }


    # Calculate Standardized Mean Difference (SMD) for overlap
    X_mat <- as.matrix(X_all[, paste0("X", 1:p_obs)])
    X_trt <- X_mat[idx_rct_trt, , drop = FALSE]
    X_ctrl <- X_mat[idx_rct_ctrl, , drop = FALSE]
    X_rwd <- X_mat[idx_rwd_ctrl, , drop = FALSE]

    for (j in 1:p_obs) {
      # Pooled SD for each comparison
      pooled_sd_trt_ctrl <- sqrt((var(X_trt[,j]) + var(X_ctrl[,j])) / 2)
      pooled_sd_trt_rwd <- sqrt((var(X_trt[,j]) + var(X_rwd[,j])) / 2)
      pooled_sd_ctrl_rwd <- sqrt((var(X_ctrl[,j]) + var(X_rwd[,j])) / 2)

      # SMD = (mean1 - mean2) / pooled_sd
      smd_trt_ctrl_all[i, j] <- (mean(X_trt[,j]) - mean(X_ctrl[,j])) / pooled_sd_trt_ctrl
      smd_trt_rwd_all[i, j] <- (mean(X_trt[,j]) - mean(X_rwd[,j])) / pooled_sd_trt_rwd
      smd_ctrl_rwd_all[i, j] <- (mean(X_ctrl[,j]) - mean(X_rwd[,j])) / pooled_sd_ctrl_rwd
    }

    # Calculate SMD for U (sc == 3 only)
    if (sc_test == 3 && !is.null(data_sample$U)) {
      U_mat <- as.matrix(data_sample$U)
      U_trt <- U_mat[idx_rct_trt, , drop = FALSE]
      U_ctrl <- U_mat[idx_rct_ctrl, , drop = FALSE]
      U_rwd <- U_mat[idx_rwd_ctrl, , drop = FALSE]

      for (j in 1:n_U) {
        pooled_sd_trt_ctrl_U <- sqrt((var(U_trt[,j]) + var(U_ctrl[,j])) / 2)
        pooled_sd_trt_rwd_U <- sqrt((var(U_trt[,j]) + var(U_rwd[,j])) / 2)
        pooled_sd_ctrl_rwd_U <- sqrt((var(U_ctrl[,j]) + var(U_rwd[,j])) / 2)

        smd_U_trt_ctrl_all[i, j] <- (mean(U_trt[,j]) - mean(U_ctrl[,j])) / pooled_sd_trt_ctrl_U
        smd_U_trt_rwd_all[i, j] <- (mean(U_trt[,j]) - mean(U_rwd[,j])) / pooled_sd_trt_rwd_U
        smd_U_ctrl_rwd_all[i, j] <- (mean(U_ctrl[,j]) - mean(U_rwd[,j])) / pooled_sd_ctrl_rwd_U
      }
    }


    # Print summary statistics for first iteration only
    if (i == 1 && c_test == cor[1]) {
      cat("\n=== Summary Statistics ===\n")
      cat("Sample size: n_rct_trt =", sum(D == 1 & Z == 1),
          ", n_rct_ctrl =", sum(D == 1 & Z == 0),
          ", n_rwd =", sum(D == 0), "\n")
      cat("True ATE:", data_sample$treat_eff, "\n")
      cat("True population mean (treatment):", data_sample$true_mean_trt_pop, "\n")
      cat("True population mean (control):", data_sample$true_mean_ctrl_pop, "\n")
      cat("RCT sigma:", data_sample$sigma_rct, "\n")
      cat("RWD sigma:", data_sample$sigma_rwd, "\n")
    }
  }

  cat("\n--- Correlation =", c_test, "---\n")
  if (sc_test == 3) {
    cat("R² using U (true confounders):\n")
    cat("  RCT Treatment (D=1, Z=1): R² =", round(mean(R2_rct_trt_all), 4), "\n")
    cat("  RCT Control (D=1, Z=0):   R² =", round(mean(R2_rct_ctrl_all), 4), "\n")
    cat("  RWD Control (D=0, Z=0):   R² =", round(mean(R2_rwd_ctrl_all), 4), "\n")
    cat("R² using X (observed covariates):\n")
    cat("  RCT Treatment (D=1, Z=1): R² =", round(mean(R2_rct_trt_X_all), 4), "\n")
    cat("  RCT Control (D=1, Z=0):   R² =", round(mean(R2_rct_ctrl_X_all), 4), "\n")
    cat("  RWD Control (D=0, Z=0):   R² =", round(mean(R2_rwd_ctrl_X_all), 4), "\n")
  } else {
    cat("RCT Treatment (D=1, Z=1): R² =", round(mean(R2_rct_trt_all), 4), "\n")
    cat("RCT Control (D=1, Z=0):   R² =", round(mean(R2_rct_ctrl_all), 4), "\n")
    cat("RWD Control (D=0, Z=0):   R² =", round(mean(R2_rwd_ctrl_all), 4), "\n")
  }

  # Print overlap measures for X
  cat("\nOverlap Measures for X (Mean Absolute SMD across simulations):\n")
  mean_abs_smd_trt_ctrl <- colMeans(abs(smd_trt_ctrl_all))
  mean_abs_smd_trt_rwd <- colMeans(abs(smd_trt_rwd_all))
  mean_abs_smd_ctrl_rwd <- colMeans(abs(smd_ctrl_rwd_all))

  cat("  RCT Trt vs RCT Ctrl: ", paste0("X", 1:p_obs, "=", round(mean_abs_smd_trt_ctrl, 3), collapse = ", "), "\n")
  cat("    Mean ASMD =", round(mean(mean_abs_smd_trt_ctrl), 4), "\n")
  cat("  RCT Trt vs RWD Ctrl: ", paste0("X", 1:p_obs, "=", round(mean_abs_smd_trt_rwd, 3), collapse = ", "), "\n")
  cat("    Mean ASMD =", round(mean(mean_abs_smd_trt_rwd), 4), "\n")
  cat("  RCT Ctrl vs RWD Ctrl:", paste0("X", 1:p_obs, "=", round(mean_abs_smd_ctrl_rwd, 3), collapse = ", "), "\n")
  cat("    Mean ASMD =", round(mean(mean_abs_smd_ctrl_rwd), 4), "\n")

  # Print overlap measures for U (sc == 3 only)
  if (sc_test == 3) {
    cat("\nOverlap Measures for U (Mean Absolute SMD across simulations):\n")
    mean_abs_smd_U_trt_ctrl <- colMeans(abs(smd_U_trt_ctrl_all), na.rm = TRUE)
    mean_abs_smd_U_trt_rwd <- colMeans(abs(smd_U_trt_rwd_all), na.rm = TRUE)
    mean_abs_smd_U_ctrl_rwd <- colMeans(abs(smd_U_ctrl_rwd_all), na.rm = TRUE)

    cat("  RCT Trt vs RCT Ctrl: ", paste0("U", 1:n_U, "=", round(mean_abs_smd_U_trt_ctrl, 3), collapse = ", "), "\n")
    cat("    Mean ASMD =", round(mean(mean_abs_smd_U_trt_ctrl), 4), "\n")
    cat("  RCT Trt vs RWD Ctrl: ", paste0("U", 1:n_U, "=", round(mean_abs_smd_U_trt_rwd, 3), collapse = ", "), "\n")
    cat("    Mean ASMD =", round(mean(mean_abs_smd_U_trt_rwd), 4), "\n")
    cat("  RCT Ctrl vs RWD Ctrl:", paste0("U", 1:n_U, "=", round(mean_abs_smd_U_ctrl_rwd, 3), collapse = ", "), "\n")
    cat("    Mean ASMD =", round(mean(mean_abs_smd_U_ctrl_rwd), 4), "\n")
  }

}
