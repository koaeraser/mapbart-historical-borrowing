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

size_option <- "small"

if (size_option == "small") {
  n1      <- 30
} else {
  n1      <- (r/(r+1)) * n
}
n3 <- n
size_suffix <- paste0("_n", n1)

n_rwd_pool <- 2 * (n1 + n3)
n_rct_pool <- 2 * (n1 + n3)

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

seed_rwd <- 43

# RWD-generation mode: TRUE = one frozen pool (files tagged _fz); FALSE = fresh RWD each iteration
rwd_frozen <- FALSE
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
  eff <- log(1.4)
}
eff_star <- log(1.0)

if (subsc == "E"){
  if (sc == 1) {
    b0_rct <- 0.20
    b0_rwd <- 0.20

    sd_Y_rct <- 1.25
    sd_Y_rwd <- 1.25

  } else if (sc == 2) {
    b0_rct <- 0.10
    b0_rwd <- 0.10

    sd_Y_rct <- 1.25
    sd_Y_rwd <- 1.25

  } else if (sc == 3) {
    b0_rct <- 0.10
    b0_rwd <- 0.10

    sd_Y_rct <- 1.25
    sd_Y_rwd <- 1.25

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
# Assigns the pool objects (X/W/y/lp/yobs/censor/label/C/entry/admin _rwd_pool) to global scope.
gen_rwd_pool <- function(seed_x, seed_y) {
  set.seed(seed_x)
  Xp <- matrix(rep(NA, n_rwd_pool * p_obs), ncol = p_obs)
  if (sc == 3) {
    Up <- matrix(rep(NA, n_rwd_pool * length(conf_cols)), ncol = length(conf_cols))
  }
  for (i in 1:n_rwd_pool) {
    #========================== Generate X ==========================
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
  Cp <- rexp(n_rwd_pool, rate = 0.09)        # random dropout (~12.5%)
  entryp <- runif(n_rwd_pool, 0, 0.5)
  adminp <- 5.5 - entryp                     # admin follow-up ~ U(5.0, 5.5) yr; horizon >= rmst tau (=5), ~12.5% admin censoring

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
    yp[i] <- exp(lp + rnorm(n = 1, sd = sd_Y_rwd))
  }

  yobsp <- rep(NA, n_rwd_pool)
  censorp <- rep(NA, n_rwd_pool)
  labelp <- rep(NA, n_rwd_pool)
  for (i in 1:n_rwd_pool) {
    yobsp[i] <- pmin(yp[i], Cp[i], adminp[i])
    censorp[i] <- as.integer((yp[i] <= Cp[i]) & (yp[i] <= adminp[i]))
    if (censorp[i] == 1) {
      labelp[i] <- "observed"
    } else {
      if (Cp[i] <= adminp[i]) labelp[i] <- "dropout" else labelp[i] <- "admin"
    }
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
  C_rwd_pool      <<- Cp
  entry_rwd_pool  <<- entryp
  admin_rwd_pool  <<- adminp
  W_rwd_pool      <<- Wp
  y_rwd_pool      <<- yp
  lp_rwd_pool     <<- lpp
  yobs_rwd_pool   <<- yobsp
  censor_rwd_pool <<- censorp
  label_rwd_pool  <<- labelp
  X_rwd_pool_by_c <<- X_rwd_pool_by_c
  if (sc == 3) U_rwd_pool <<- Up
}

# Frozen mode: build the pool once, up front (per-iteration mode rebuilds it inside the loop).
if (rwd_frozen) gen_rwd_pool(seed_rwd, seed_rwd + 1000)

for (ci in seq_along(cor)) {
  c <- cor[ci]
  if (rwd_frozen) X_rwd_pool <- X_rwd_pool_by_c[[ci]]

  # Per-iteration output path -- SINGLE definition, used for both the save below
  # and the iter-1 skip check, so the two can never drift.
  iter_file <- function(it)
    paste0(mainDir, "/mapbart-sim-survival-realcov/data_v3/data_p", p_obs, size_suffix,
           "_sc", sc, if (!is.null(subsc)) subsc else "",
           if (sc == 1 | sc == 2) "" else paste0("_cor", c),
           "_", hypo, if (rwd_frozen) "_fz" else "", "_", it, ".RData")

  #===========================================================================
  #============================ Generate RCT_pool ============================
  #===========================================================================
  for (iter in 1:niter){

    # Skip regenerating iteration 1 when its data file already exists (e.g. the
    # single-iteration file the ess calibration was run on).  Each iteration
    # fully resets the RNG via set.seed(seed[ee]) with no carryover, and seed[1]
    # is independent of niter, so the existing _1 is exactly what this loop would
    # regenerate -- we only advance ee to keep iterations 2..niter bit-identical.
    if (iter == 1 && file.exists(iter_file(1))) {
      cat(sprintf("data_gen: iter 1 exists (%s) -- skipping regeneration\n",
                  basename(iter_file(1))))
      ee <- ee + 1
      next
    }

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
      # Calibrated intercept
      target_sel <- n1 / (n1 + n3)
      thres <- uniroot(function(cc) mean(plogis((lp - cc) / sel_temp)) - target_sel,
                       interval = range(lp) + c(-10, 10))$root
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }
    if (sc == 3) {
      if (length(beta_D) > 1) {
        lp <- U[, 1:length(conf_cols)] %*% beta_D
      } else {
        lp <- beta_D * rowSums(U[, 1:length(conf_cols)])
      }
      # Calibrated intercept
      target_sel <- n1 / (n1 + n3)
      thres <- uniroot(function(cc) mean(plogis((lp - cc) / sel_temp)) - target_sel,
                       interval = range(lp) + c(-10, 10))$root
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }
    
    C_rct_gen <- rexp(n_rct_pool, rate = 0.09)     # random dropout (~12.5%)
    entry_rct_gen <- runif(n_rct_pool, 0, 0.5)
    admin_rct_gen <- 5.5 - entry_rct_gen           # admin follow-up ~ U(5.0, 5.5) yr; horizon >= rmst tau (=5), ~12.5% admin censoring
    
    C <- c(C_rct_gen, C_rwd_pool)
    entry <- c(entry_rct_gen, entry_rwd_pool)
    admin <- c(admin_rct_gen, admin_rwd_pool)
    
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
    C_rct_sel <- C[idx_rct]
    entry_rct_sel <- entry[idx_rct]
    admin_rct_sel <- admin[idx_rct]
    if (sc == 3) {
      W_rct_sel <- matrix(U[idx_rct, ], ncol = length(conf_cols))
    } else {
      W_rct_sel <- X_rct_sel[, conf_cols]
    }
    
    #========================== Generate Z ==========================
    Z_rct <- rep(1, n_rct)
    
    #========================== Generate Y ==========================
    y_rct <- rep(NA, n_rct)
    yobs_rct <- rep(NA, n_rct)
    label_rct <- rep(NA, n_rct)
    lp_rct <- rep(NA, n_rct)
    
    true_median_trt <- rep(NA, n_rct)
    true_median_ctrl <- rep(NA, n_rct)
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
        y_rct[i] <- exp(lp + modifier + eff_i + rnorm(n = 1, sd = sd_Y_rct))
      } else {
        lp_rct[i] <- lp + modifier
        y_rct[i] <- exp(lp + modifier + rnorm(n = 1, sd = sd_Y_rct))
      }
      
      eff_i_rct[i] <- eff_i
      true_median_trt[i] <- exp(lp + modifier + eff_i)
      true_median_ctrl[i] <- exp(lp + modifier)
    }
    
    censor_rct <- rep(NA, n_rct)
    for (i in 1:n_rct) {
      yobs_rct[i] <- pmin(y_rct[i], C_rct_sel[i], admin_rct_sel[i])
      censor_rct[i] <- as.integer((y_rct[i] <= C_rct_sel[i]) & (y_rct[i] <= admin_rct_sel[i]))
      
      if (censor_rct[i] == 1) {
        label_rct[i] <- "observed"
      } else {
        if (C_rct_sel[i] <= admin_rct_sel[i]) {
          label_rct[i] <- "dropout"
        } else {
          label_rct[i] <- "admin"
        }
      }
    }
    
    S_mix_fun <- function(t, mu_vec, sigma) {
      mean(1 - pnorm((log(t) - mu_vec) / sigma))
    }
    
    mu_trt_vec <- log(true_median_trt)
    f <- function(t) S_mix_fun(t, mu_trt_vec, sd_Y_rct) - 0.5
    median_trt_pop <- uniroot(f,
                              interval = c(exp(min(mu_trt_vec) - 3*sd_Y_rct),
                                           exp(max(mu_trt_vec) + 3*sd_Y_rct)))$root
    
    mu_ctrl_vec <- log(true_median_ctrl)
    f <- function(t) S_mix_fun(t, mu_ctrl_vec, sd_Y_rct) - 0.5
    median_ctrl_pop <- uniroot(f,
                               interval = c(exp(min(mu_ctrl_vec) - 3*sd_Y_rct),
                                            exp(max(mu_ctrl_vec) + 3*sd_Y_rct)))$root
    
    rmst_tau_dgp <- 5
    rmst_pop_fun <- function(mu_vec, sigma, tau, K = 200) {
      grid <- seq(1e-3, tau, length.out = K)
      mean(vapply(grid, function(t) S_mix_fun(t, mu_vec, sigma), numeric(1))) * tau
    }
    rmst_trt_pop  <- rmst_pop_fun(mu_trt_vec,  sd_Y_rct, rmst_tau_dgp)
    rmst_ctrl_pop <- rmst_pop_fun(mu_ctrl_vec, sd_Y_rct, rmst_tau_dgp)
    
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
    yobs_rwd <- yobs_rwd_pool[pool_indices]
    censor_rwd <- censor_rwd_pool[pool_indices]
    label_rwd <- label_rwd_pool[pool_indices]
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
        "eff_true =", log(median_trt_pop / median_ctrl_pop),"\n")
    # Output
    out <- list(X = rbind(X_rct_df, X_rwd_df),
                U = U_combined,
                y = c(yobs_rct, yobs_rwd),
                y_true = c(y_rct, y_rwd),
                delta = c(censor_rct, censor_rwd),
                label = c(label_rct, label_rwd),
                treat_eff = eff,
                treat_eff_true = log(median_trt_pop / median_ctrl_pop),
                treat_eff_star = eff_star,
                gamma_eff = gamma,  # modifier coefficient
                eta_eff = if (exists("eta")) eta else 0,  # HTE coefficient for treatment effect
                sigma_rct = sd_Y_rct,
                sigma_rwd = sd_Y_rwd,
                true_median_trt = true_median_trt,
                true_median_ctrl = true_median_ctrl,
                true_median_trt_pop = median_trt_pop,
                true_median_ctrl_pop = median_ctrl_pop,
                true_rmst_tau = rmst_tau_dgp,
                true_rmst_trt_pop = rmst_trt_pop,
                true_rmst_ctrl_pop = rmst_ctrl_pop,
                eff_i = eff_i_rct,
                lp = c(lp_rct, lp_rwd),
                n_rct = n_rct,
                n_rwd = n_rwd
    )
    
    # Save with scenario identifiers (path from iter_file so it matches the skip check)
    saveRDS(out, file = iter_file(iter))
  }
}

## Testing code - update sc and subsc as needed
sc_test <- sc
subsc_test <- subsc

cat("\n========================================\n")
cat("Testing R² for sc =", sc_test, ", subsc =", ifelse(is.null(subsc_test), "NULL", subsc_test), "\n")
cat("========================================\n")

for (c_test in cor) {
  
  rate1 <- rate2 <- rate3 <- 0
  rate_admin1 <- rate_admin2 <- rate_admin3 <- 0
  rate_drop1 <- rate_drop2 <- rate_drop3 <- 0
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


  # Median survival times across iterations
  median_surv_rct_trt_all <- numeric(niter)
  median_surv_rct_ctrl_all <- numeric(niter)
  median_surv_rwd_ctrl_all <- numeric(niter)
  pop_median_trt_all <- numeric(niter)
  pop_median_ctrl_all <- numeric(niter)
  
  for (i in 1:niter){
    if (sc_test == 1 | sc_test == 2) {
      if (is.null(subsc_test)) {
        filepath <- paste0(mainDir,"/mapbart-sim-survival-realcov/data_v3/data_p",p_obs,size_suffix,"_sc",sc_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      } else {
        filepath <- paste0(mainDir,"/mapbart-sim-survival-realcov/data_v3/data_p",p_obs,size_suffix,"_sc",sc_test,subsc_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      }
    } else {
      if (is.null(subsc_test)) {
        filepath <- paste0(mainDir,"/mapbart-sim-survival-realcov/data_v3/data_p",p_obs,size_suffix,"_sc",sc_test,"_cor",c_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      } else {
        filepath <- paste0(mainDir,"/mapbart-sim-survival-realcov/data_v3/data_p",p_obs,size_suffix,"_sc",sc_test,subsc_test,"_cor",c_test,"_",hypo, if (rwd_frozen) "_fz" else "", "_",i,".RData")
      }
    }

    # Read data
    data_sample <- readRDS(filepath)
    delta <- data_sample$delta
    label <- data_sample$label
    X_all <- data_sample$X
    D <- X_all[,"D"]
    Z <- X_all[,"Z"]
    
    # Calculate event rates
    tab1 <- table(label[D == 1 & Z == 1])
    tab2 <- table(label[D == 1 & Z == 0])
    tab3 <- table(label[D == 0 & Z == 0])
    
    rate1 = rate1 + ifelse("observed" %in% names(tab1), tab1["observed"], 0) / length(label[D == 1 & Z == 1])
    rate2 = rate2 + ifelse("observed" %in% names(tab2), tab2["observed"], 0) / length(label[D == 1 & Z == 0])
    rate3 = rate3 + ifelse("observed" %in% names(tab3), tab3["observed"], 0) / length(label[D == 0 & Z == 0])
    
    rate_admin1 = rate_admin1 + ifelse("admin" %in% names(tab1), tab1["admin"], 0) / length(label[D == 1 & Z == 1])
    rate_admin2 = rate_admin2 + ifelse("admin" %in% names(tab2), tab2["admin"], 0) / length(label[D == 1 & Z == 0])
    rate_admin3 = rate_admin3 + ifelse("admin" %in% names(tab3), tab3["admin"], 0) / length(label[D == 0 & Z == 0])
    
    rate_drop1 = rate_drop1 + ifelse("dropout" %in% names(tab1), tab1["dropout"], 0) / length(label[D == 1 & Z == 1])
    rate_drop2 = rate_drop2 + ifelse("dropout" %in% names(tab2), tab2["dropout"], 0) / length(label[D == 1 & Z == 0])
    rate_drop3 = rate_drop3 + ifelse("dropout" %in% names(tab3), tab3["dropout"], 0) / length(label[D == 0 & Z == 0])
    
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
    
    # Calculate median survival times for this iteration
    Y_iter <- if (!is.null(data_sample$y_true)) data_sample$y_true else data_sample$y
    median_surv_rct_trt_all[i] <- median(Y_iter[idx_rct_trt])
    median_surv_rct_ctrl_all[i] <- median(Y_iter[idx_rct_ctrl])
    median_surv_rwd_ctrl_all[i] <- median(Y_iter[idx_rwd_ctrl])
    pop_median_trt_all[i] <- data_sample$true_median_trt_pop
    pop_median_ctrl_all[i] <- data_sample$true_median_ctrl_pop

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

      Y_all <- if (!is.null(data_sample$y_true)) data_sample$y_true else data_sample$y

      eff_i_X <- eff

      lp_X_rct_trt <- lp_X_rct[idx_rct_trt] + gamma * rowSums(X_mat[idx_rct_trt, ]) + eff_i_X
      Y_rct_trt <- log(Y_all[idx_rct_trt])
      R2_rct_trt_X_all[i] <- cor(Y_rct_trt, lp_X_rct_trt)^2

      lp_X_rct_ctrl <- lp_X_rct[idx_rct_ctrl] + gamma * rowSums(X_mat[idx_rct_ctrl, ])
      Y_rct_ctrl <- log(Y_all[idx_rct_ctrl])
      R2_rct_ctrl_X_all[i] <- cor(Y_rct_ctrl, lp_X_rct_ctrl)^2

      lp_X_rwd_ctrl <- lp_X_rwd[idx_rwd_ctrl]
      Y_rwd_ctrl <- log(Y_all[idx_rwd_ctrl])
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

  }

  cat("\n--- Correlation =", c_test, "---\n")
  
  # Print mean of median survival times averaged over all iterations
  cat("Mean of Median Survival Time (true, averaged over", niter, "iterations):\n")
  cat("  RCT Trt  (D=1, Z=1):", round(mean(median_surv_rct_trt_all), 4), "\n")
  cat("  RCT Ctrl (D=1, Z=0):", round(mean(median_surv_rct_ctrl_all), 4), "\n")
  cat("  RWD Ctrl (D=0, Z=0):", round(mean(median_surv_rwd_ctrl_all), 4), "\n")
  cat("Mean of Population Median (true, averaged over", niter, "iterations):\n")
  cat("  Trt Pop Median:", round(mean(pop_median_trt_all), 4), "\n")
  cat("  Ctrl Pop Median:", round(mean(pop_median_ctrl_all), 4), "\n")
  
  cat("Event Rates: RCT_trt =", round(rate1 / niter, 4), ", RCT_ctrl =", round(rate2 / niter, 4), ", RWD =", round(rate3 / niter, 4), "\n")
  cat("Admin Censoring: RCT_trt =", round(rate_admin1 / niter, 4), ", RCT_ctrl =", round(rate_admin2 / niter, 4), ", RWD =", round(rate_admin3 / niter, 4), "\n")
  cat("Dropout Censoring: RCT_trt =", round(rate_drop1 / niter, 4), ", RCT_ctrl =", round(rate_drop2 / niter, 4), ", RWD =", round(rate_drop3 / niter, 4), "\n")
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
