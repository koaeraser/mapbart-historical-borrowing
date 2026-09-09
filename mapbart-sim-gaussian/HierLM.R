rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop/"
library(rstan)

data_folder <- "data_v2"  # Options: "data" or "data_v2"

sc <- 3
subsc <- "E"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-gaussian", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-gaussian", data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# prior_vals <- c(0.1, 0.25, 0.5, 1.0, 5.0, 10.0)
prior_vals <- 0.25

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
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, subsc, "_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, "_cor", cor, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, subsc, "_cor", cor, "_threshold.RData")
    }
  }

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}

# Stan sampling parameters
n_chains <- 1
n_iter <- 1600
n_warmup <- 100
n_cores <- 3

# Options for mu_alpha and mu_beta (hierarchical hyperprior means)
estimate_mu_alpha <- TRUE
estimate_mu_beta <- TRUE
mu_alpha_fixed <- 0   # used when estimate_mu_alpha = FALSE
mu_beta_fixed <- 0    # used when estimate_mu_beta = FALSE

# Define Stan model for Hierarchical Bayesian Linear Regression (for control group)
lm_hierarchical_code <- "
data {
  int<lower=0> N;           // number of observations
  int<lower=0> P;           // number of predictors
  matrix[N, P] X;           // predictor matrix
  vector[N] y;              // response variable
  int<lower=0,upper=1> D[N]; // data source indicator (0=RWD, 1=RCT)
  vector<lower=0>[2] sigma_scale; // scale parameters for IG prior on sigma2 (1=RWD, 2=RCT)
  int<lower=0,upper=1> estimate_mu_alpha; // 1 = estimate mu_alpha, 0 = use fixed value
  int<lower=0,upper=1> estimate_mu_beta;  // 1 = estimate mu_beta, 0 = use fixed value
  real mu_alpha_fixed; // fixed value for mu_alpha (used when estimate_mu_alpha = 0)
  real mu_beta_fixed;  // fixed value for mu_beta (used when estimate_mu_beta = 0)
  real<lower=0> prior;      // prior scale for tau2 ~ inv_gamma(3/2, 3*prior/2)
}

parameters {
  vector[2] alpha;          // intercept for each group (1=RWD, 2=RCT)
  real<lower=0> tau_a2;     // variance for hierarchical prior on alpha
  matrix[2, P] beta;        // regression coefficients for each group (row 1=RWD, row 2=RCT)
  vector<lower=0>[P] tau2;  // predictor-specific variance for hierarchical prior on beta
  vector<lower=0>[2] sigma2; // residual variance for each group (1=RWD, 2=RCT)
  real mu_alpha;            // hyperprior mean for alpha (estimated or fixed)
  real mu_beta;             // hyperprior mean for beta (estimated or fixed)
}

model {
  vector[N] mu;
  vector[2] sigma;
  real tau_a;
  vector[P] tau;

  // Hyperprior means (estimated or pinned to a fixed value)
  if (estimate_mu_alpha == 1) {
    mu_alpha ~ normal(0, 10);
  } else {
    mu_alpha ~ normal(mu_alpha_fixed, 0.0001);
  }
  if (estimate_mu_beta == 1) {
    mu_beta ~ normal(0, 10);
  } else {
    mu_beta ~ normal(mu_beta_fixed, 0.0001);
  }

  // Hierarchical prior on intercepts
  tau_a2 ~ inv_gamma(3.0/2.0, 3.0*prior/2.0);
  tau_a = sqrt(tau_a2);
  for (g in 1:2) {
    alpha[g] ~ normal(mu_alpha, tau_a);
  }

  // Hierarchical prior on coefficients (predictor-specific variance)
  for (j in 1:P) {
    tau2[j] ~ inv_gamma(3.0/2.0, 3.0*prior/2.0);
    tau[j] = sqrt(tau2[j]);
  }
  for (g in 1:2) {
    for (j in 1:P) {
      beta[g, j] ~ normal(mu_beta, tau[j]);
    }
  }

  // Priors for group-specific residual variances (BART-style)
  for (g in 1:2) {
    sigma2[g] ~ inv_gamma(3.0/2.0, 3.0*sigma_scale[g]/2.0);
  }

  // Derived parameters
  for (g in 1:2) {
    sigma[g] = sqrt(sigma2[g]);
  }

  // Likelihood - use appropriate coefficients and variance based on group
  for (n in 1:N) {
    int group_idx = D[n] + 1;  // Convert 0/1 to 1/2 for indexing
    mu[n] = alpha[group_idx] + X[n] * beta[group_idx]';
    y[n] ~ normal(mu[n], sigma[group_idx]);
  }
}
"

# Compile Stan model for control (hierarchical)
lm_hierarchical_model <- stan_model(model_code = lm_hierarchical_code)

# Define simple Stan model for treatment (like LMv1)
lm_simple_code <- "
data {
  int<lower=0> N;           // number of observations
  int<lower=0> P;           // number of predictors
  matrix[N, P] X;           // predictor matrix
  vector[N] y;              // response variable
  real<lower=0> sigma_scale; // scale parameter for IG prior on sigma2
  real mu_alpha;            // prior mean for intercept
  real<lower=0> lambda_alpha; // BART-style prior SD for intercept
  real<lower=0> lambda_beta;  // BART-style prior SD for coefficients
}

parameters {
  real alpha;               // intercept
  vector[P] beta;           // regression coefficients
  real<lower=0> sigma2;     // variance parameter
}

model {
  vector[N] mu;
  real sigma;

  // Priors (BART-style, matching the AFT models)
  alpha ~ normal(mu_alpha, lambda_alpha);
  beta ~ normal(0, lambda_beta);
  sigma2 ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0); // BART-style prior

  // Derived parameter
  sigma = sqrt(sigma2);

  // Linear predictor
  mu = alpha + X * beta;

  // Likelihood for linear model
  y ~ normal(mu, sigma);
}
"

# Compile simple Stan model for treatment
lm_simple_model <- stan_model(model_code = lm_simple_code)

for (prior in prior_vals) {

  cat("\n========== Running with prior =", prior, "==========\n")

  set.seed(6)
  seed <- sample(1:10000, niter*4, replace = F)
  ee <- 1

  res <- data.frame(
    dist = "normal",
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
  res$pehe_subj <- NA
  res$bias_subj <- NA
  res$bias.trt.sigma <- NA
  res$sd.trt.sigma <- NA
  res$w2distance.trt.sigma <- NA
  res$bias.ctrl.sigma <- NA
  res$sd.ctrl.sigma <- NA
  res$w2distance.ctrl.sigma <- NA
  res$bias.trt.median.pop <- NA
  res$w2distance.trt.median.pop <- NA
  res$bias.ctrl.median.pop <- NA
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
    this_seed <- seed[ee]   # one seed per iteration; both arms share it
    set.seed(this_seed)
    ee <- ee + 1


    data = list()

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
    data$X <- rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
                    data_rwd_full$X[rwd_idx, paste0("X",1:10)])
    data$D <- c(data_rct_full$X[rct_idx, "D"],
                data_rwd_full$X[rwd_idx, "D"])
    data$Z <- c(data_rct_full$X[rct_idx, "Z"],
                data_rwd_full$X[rwd_idx, "Z"])
    data$Y <- c(data_rct_full$y[rct_idx],
                data_rwd_full$y[rwd_idx])

    # Also need filename variable for reading true values later
    filename <- filename_rct

    #------------------------------------------------
    #----------- RCT + Historical Control -----------
    #------------------------------------------------
    # Center X at the RCT/target (D==1) covariate means (common reference, matching AFT models)
    x_means <- colMeans(as.matrix(data$X[data$D == 1, ]))
    # Pool all control data (Z=0) from both RCT (D=1) and RWD (D=0)
    x_ctrl <- as.matrix(data$X[data$Z == 0, ])
    x_ctrl <- scale(x_ctrl, center = x_means, scale = FALSE)
    y_ctrl <- data$Y[data$Z == 0]
    d_ctrl <- data$D[data$Z == 0]  # Data source indicator

    # Estimate sigma separately for RWD and RCT groups (BART approach)
    nu <- 3
    sigquant <- 0.90
    qchi <- qchisq(1 - sigquant, nu)

    # For RWD group (D=0)
    x_ctrl_rwd <- x_ctrl[d_ctrl == 0, ]
    y_ctrl_rwd <- y_ctrl[d_ctrl == 0]
    if (length(y_ctrl_rwd) > 0) {
      if (ncol(x_ctrl_rwd) < length(y_ctrl_rwd)) {
        lmf_rwd <- lm(y_ctrl_rwd ~ x_ctrl_rwd)
        sigest_rwd <- summary(lmf_rwd)$sigma
      } else {
        sigest_rwd <- sd(y_ctrl_rwd)
      }
    } else {
      sigest_rwd <- sd(y_ctrl)  # Fallback if no RWD data
    }
    lambda_rwd <- (sigest_rwd^2 * qchi) / nu

    # For RCT group (D=1)
    x_ctrl_rct <- x_ctrl[d_ctrl == 1, ]
    y_ctrl_rct <- y_ctrl[d_ctrl == 1]
    if (length(y_ctrl_rct) > 0) {
      if (ncol(x_ctrl_rct) < length(y_ctrl_rct)) {
        lmf_rct <- lm(y_ctrl_rct ~ x_ctrl_rct)
        sigest_rct <- summary(lmf_rct)$sigma
      } else {
        sigest_rct <- sd(y_ctrl_rct)
      }
    } else {
      sigest_rct <- sd(y_ctrl)  # Fallback if no RCT control data
    }
    lambda_rct <- (sigest_rct^2 * qchi) / nu

    # Prepare data for Stan with hierarchical structure
    stan_data_ctrl <- list(
      N = length(y_ctrl),
      P = ncol(x_ctrl),
      X = x_ctrl,
      y = y_ctrl,
      D = d_ctrl,  # Pass data source indicator
      sigma_scale = c(lambda_rwd, lambda_rct),  # Pass group-specific scale parameters
      estimate_mu_alpha = as.integer(estimate_mu_alpha),
      estimate_mu_beta = as.integer(estimate_mu_beta),
      mu_alpha_fixed = mu_alpha_fixed,
      mu_beta_fixed = mu_beta_fixed,
      prior = prior  # Prior scale for tau2
    )

    # Fit hierarchical linear model for control group
    fit_ctrl <- sampling(lm_hierarchical_model,
                         data = stan_data_ctrl,
                         chains = n_chains,
                         iter = n_iter,
                         warmup = n_warmup,
                         cores = n_cores,
                         seed = this_seed,
                         refresh = 0)

    # Extract posterior samples
    samples_ctrl <- rstan::extract(fit_ctrl)

    # Predict for target population (D=1) using RCT coefficients
    x_test <- as.matrix(data$X[data$D == 1, ])
    x_test <- scale(x_test, center = x_means, scale = FALSE)

    n_samples <- length(samples_ctrl$alpha[,1])
    n_test <- nrow(x_test)

    # Generate predictions using RCT group coefficients (group_idx = 2)
    mu_pred_ctrl <- matrix(NA, n_samples, n_test)
    for (i in 1:n_samples) {
      mu_pred_ctrl[i, ] <- samples_ctrl$alpha[i, 2] + x_test %*% samples_ctrl$beta[i, 2, ]
    }

    #------------------------------------------------
    #---------------- RCT Treatment -----------------
    #------------------------------------------------
    x_trt <- as.matrix(data$X[data$Z == 1, ])
    x_trt <- scale(x_trt, center = x_means, scale = FALSE)
    y_trt <- data$Y[data$Z == 1]

    # Estimate sigma for treatment group (BART approach)
    if (ncol(x_trt) < length(y_trt)) {
      lmf_trt <- lm(y_trt ~ x_trt)
      sigest_trt <- summary(lmf_trt)$sigma
    } else {
      sigest_trt <- sd(y_trt)
    }
    lambda_trt <- (sigest_trt^2 * qchi) / nu

    # BART-style coefficient prior scales (range budget, alpha-weighted) -- matches AFT
    k <- 2.0
    total_var_trt <- (max(y_trt) - min(y_trt))^2 / (2 * k)^2
    w_alpha <- 2.0
    w_beta <- 1.0
    total_weight_trt <- w_alpha + ncol(x_trt) * w_beta
    lambda_alpha_trt <- sqrt(total_var_trt * w_alpha / total_weight_trt)
    lambda_beta_trt  <- sqrt(total_var_trt * w_beta  / total_weight_trt)

    # Prepare data for Stan (simple model)
    stan_data_trt <- list(
      N = length(y_trt),
      P = ncol(x_trt),
      X = x_trt,
      y = y_trt,
      sigma_scale = lambda_trt,  # Single scale parameter
      mu_alpha = 0,
      lambda_alpha = lambda_alpha_trt,
      lambda_beta = lambda_beta_trt
    )

    # Fit simple linear model for treatment group (like LMv1)
    fit_trt <- sampling(lm_simple_model,
                        data = stan_data_trt,
                        chains = n_chains,
                        iter = n_iter,
                        warmup = n_warmup,
                        cores = n_cores,
                        seed = this_seed,
                        refresh = 0)

    # Extract posterior samples
    samples_trt <- rstan::extract(fit_trt)

    # Predict for target population (D=1)
    n_samples_trt <- length(samples_trt$alpha)

    # Generate predictions
    mu_pred_trt <- matrix(NA, n_samples_trt, n_test)
    for (i in 1:n_samples_trt) {
      mu_pred_trt[i, ] <- samples_trt$alpha[i] + x_test %*% samples_trt$beta[i, ]
    }

    #-------------------------------------
    #---------- Calculate ATE ------------
    #-------------------------------------
    post_samples <- rowMeans(mu_pred_trt - mu_pred_ctrl)
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

    res[which(res$c==c & res$iteration==iter),
        "bias"] <- bias

    res[which(res$c==c & res$iteration==iter),
        "sd"] <- SD

    res[which(res$c==c & res$iteration==iter),
        "rmse"] <- rmse

    res[which(res$c==c & res$iteration==iter),
        "w1distance"] <- w1distance

    res[which(res$c==c & res$iteration==iter),
        "w2distance"] <- w2distance

    res[which(res$c==c & res$iteration==iter),
        "ci"] <- ci

    res[which(res$c==c & res$iteration==iter),
        "coverage"] <- coverage

    res[which(res$c==c & res$iteration==iter),
        "tp_calibrated"] <- tp_calibrated

    res[which(res$c==c & res$iteration==iter),
        "tp"] <- tp

    res[which(res$c==c & res$iteration==iter),
        "fp"] <- fp

    #------------------------------------------
    #------ Calculate subject-wise ------------
    #------------------------------------------
    post_samples_i <- mu_pred_trt - mu_pred_ctrl   # n_draws x n_subjects
    eff_i <- readRDS(filename)$eff_i

    # Posterior estimates per subject
    delta_hat_i <- apply(post_samples_i, 2, mean)    # n_subjects

    # Within replicate metrics
    bias_subj <- mean(delta_hat_i - eff_i)
    pehe_subj <- sqrt(mean((delta_hat_i - eff_i)^2))

    res[which(res$c==c & res$iteration==iter), "bias_subj"] <- bias_subj
    res[which(res$c==c & res$iteration==iter), "pehe_subj"] <- pehe_subj

    #-------------------------------------
    #------ Variance estimation ----------
    #-------------------------------------
    true_sigma_trt <- true_sigma_ctrl <- readRDS(filename)$sigma_rct

    # Treatment group (simple model)
    sigma_samples_trt <- sqrt(samples_trt$sigma2)
    sigma_est_trt <- mean(sigma_samples_trt)
    sd_trt_sigma <- sd(sigma_samples_trt)
    w2distance_trt_sigma <- sqrt(mean((sigma_samples_trt - true_sigma_trt)^2, na.rm = TRUE))

    res[which(res$c==c & res$iteration==iter),
        "bias.trt.sigma"] <- sigma_est_trt - true_sigma_trt
    res[which(res$c==c & res$iteration==iter),
        "sd.trt.sigma"] <- sd_trt_sigma
    res[which(res$c==c & res$iteration==iter),
        "w2distance.trt.sigma"] <- w2distance_trt_sigma

    # Control group (uses RCT variance, group_idx=2)
    sigma_samples_ctrl <- sqrt(samples_ctrl$sigma2[, 2])
    sigma_est_ctrl <- mean(sigma_samples_ctrl)
    sd_ctrl_sigma <- sd(sigma_samples_ctrl)
    w2distance_ctrl_sigma <- sqrt(mean((sigma_samples_ctrl - true_sigma_trt)^2, na.rm = TRUE))

    res[which(res$c==c & res$iteration==iter),
        "bias.ctrl.sigma"] <- sigma_est_ctrl - true_sigma_ctrl
    res[which(res$c==c & res$iteration==iter),
        "sd.ctrl.sigma"] <- sd_ctrl_sigma
    res[which(res$c==c & res$iteration==iter),
        "w2distance.ctrl.sigma"] <- w2distance_ctrl_sigma

    #--------------------------------------------------
    #------ Population mean estimation ----------------
    #--------------------------------------------------
    true_mean_trt_pop <- readRDS(filename)$true_mean_trt_pop
    true_mean_ctrl_pop <- readRDS(filename)$true_mean_ctrl_pop

    # Calculate population means
    mean_trt_pop_hat <- mean(rowMeans(mu_pred_trt))
    mean_ctrl_pop_hat <- mean(rowMeans(mu_pred_ctrl))

    # Calculate metrics
    mu_pred_trt_samples <- rowMeans(mu_pred_trt)
    mu_pred_ctrl_samples <- rowMeans(mu_pred_ctrl)
    bias_trt_median_pop <- mean_trt_pop_hat - true_mean_trt_pop
    w2distance_trt_median_pop <- sqrt(mean((mu_pred_trt_samples - true_mean_trt_pop)^2))
    bias_ctrl_median_pop <- mean_ctrl_pop_hat - true_mean_ctrl_pop
    w2distance_ctrl_median_pop <- sqrt(mean((mu_pred_ctrl_samples - true_mean_ctrl_pop)^2))

    res[which(res$c==c & res$iteration==iter),
        "bias.trt.median.pop"] <- bias_trt_median_pop
    res[which(res$c==c & res$iteration==iter),
        "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
    res[which(res$c==c & res$iteration==iter),
        "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
    res[which(res$c==c & res$iteration==iter),
        "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop

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
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, "_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, subsc, "_threshold.RData")
      }
    } else {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, "_cor", c, "_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian/res/HierLM_p", p_obs, "_sc", sc, subsc, "_cor", c, "_threshold.RData")
      }
    }
    saveRDS(calibrated_threshold_grid, file = threshold_file)
  }

  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  # Save results with prior value in filename
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/HierLM_p",p_obs,"_sc",sc,"_prior",prior,"_",hypo,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/HierLM_p",p_obs,"_sc",sc,subsc,"_prior",prior,"_",hypo,".RData"))
    }
  } else {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/HierLM_p",p_obs,"_sc",sc,"_cor",c,"_prior",prior,"_",hypo,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian/res/HierLM_p",p_obs,"_sc",sc,subsc,"_cor",c,"_prior",prior,"_",hypo,".RData"))
    }
  }

} # End of prior loop
