rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop/"
library(rstan)
library(survival)
library(dplyr)

source(paste0(mainDir, "/mapbart-sim-survival-realcov/rmst_helpers.R"))
rmst_tau <- 5   # RMST restriction horizon (admin-censoring horizon); additional output
rmst_control_sigma <- "own"  # control-arm RMST sigma: "trt"|"own"|"adaptive" (own = control's own sigma-hat; avoids inflated trt sigma at small n)

data_folder <- "data_v3"  # Options: "data" or "data_v3"
size_option <- "default"   # "default" (RCT trt 200) or "small" (RCT trt 30)
size_suffix <- if (size_option == "small") "_n30" else "_n200"

# Power-prior downweighting of the RWD (external control) likelihood.
# rwd_w = fractional weight on each RWD observation; effective RWD sample size
# ~ rwd_w * n_RWD.  To match a MAP-BART borrowing target N: rwd_w = N / n_RWD.
# rwd_w = 1 reproduces full (unweighted) borrowing.  Output files are tagged _w<rwd_w>.
#
# rwd_w_vals is now DERIVED from the ess_plot.R targets tables (the N column),
# so AFTv2 is compared at the same effective RWD sample sizes MAP-BART borrows.
# Computed below (rwd_w_block), AFTER the scenario (sc/subsc/cor/size) is set,
# as a per-scenario union of the target Ns over all control configs.
n_rwd_nominal <- 300  # approximate RWD control count (n3 from data_gen_p10_v3.R)

sc <- 1
subsc <- "E"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-survival-realcov", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-survival-realcov", data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

threshold <- 0.95

# Strength of unmeasured confounding
# if (sc == 1 | sc == 2){                # (old)
if (sc == 1 | sc == 2){
  cor <- 1
}
if (sc == 3){
  cor <- 0.7
}

# rwd_w_block: derive rwd_w_vals from the ess_plot.R targets tables.  For THIS
# scenario (size + sc/subsc/cor), union the N column over every control config's
# _tab.RData, then rwd_w = N / n_rwd_nominal (plus 1 = full borrowing).  The
# per-config tables share the same N_target, so the union is the set of target
# Ns MAP-BART borrows to in this scenario.  Requires ess_plot.R to have run.
{
  .cal_dir    <- file.path(mainDir, "mapbart-sim-survival-realcov", "ess_local", "res")
  .size_tag   <- sub("^_", "", size_suffix)                     # "n30" / "n200"
  .sc_cal_tag <- if (sc == 1 | sc == 2)
      paste0("sc", sc, if (!is.null(subsc)) subsc else "") else
      paste0("sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor)
  .tgt_files  <- list.files(.cal_dir, full.names = TRUE,
                            pattern = paste0("^ess_calibration_.*_", .size_tag,
                                             "_nt[0-9]+_k[0-9.]+_", .sc_cal_tag,
                                             "_tab\\.RData$"))
  if (length(.tgt_files)) {
    .Ns <- unique(unlist(lapply(.tgt_files, function(f) as.integer(readRDS(f)$N))))
    rwd_w_vals <- sort(unique(c(1, round(.Ns / n_rwd_nominal, 4))), decreasing = TRUE)
    cat(sprintf("AFTv2.R: rwd_w_vals from %d targets table(s) [%s/%s], target N = {%s} -> rwd_w = %s\n",
                length(.tgt_files), .size_tag, .sc_cal_tag,
                paste(sort(.Ns, decreasing = TRUE), collapse = ", "),
                paste(rwd_w_vals, collapse = ", ")))
  } else {
    rwd_w_vals <- c(1)
    warning(sprintf(paste0("AFTv2.R: no scenario/size-matched targets table ",
                    "(ess_calibration_*_%s_nt<ntree>_k<k>_%s_tab.RData in %s) -- ",
                    "using rwd_w_vals = 1 (full borrowing only).  Run ess_plot.R first."),
                    .size_tag, .sc_cal_tag, .cal_dir))
  }
}

# Infer niter from the number of simulated data files matching this scenario
data_prefix <- if (sc == 1 | sc == 2) {
  paste0("data_p", p_obs, size_suffix, "_sc", sc, if (!is.null(subsc)) subsc else "", "_", hypo, "_")
} else {
  paste0("data_p", p_obs, size_suffix, "_sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor, "_", hypo, "_")
}
niter <- length(Sys.glob(file.path(mainDir, "mapbart-sim-survival-realcov", data_folder, paste0(data_prefix, "*.RData"))))
if (niter == 0) stop("No simulated data files found for this scenario")

# Stan sampling parameters
n_chains <- 1
n_iter <- 1600
n_warmup <- 100
n_cores <- 3

aft_stan_code <- "
data {
  int<lower=0> N;           // number of observations
  int<lower=0> P;           // number of predictors
  matrix[N, P] X;           // predictor matrix
  vector[N] y;              // log survival times
  int<lower=0,upper=1> event[N]; // event indicator
  vector<lower=0>[N] wt;    // per-obs likelihood weight (power-prior downweighting of RWD)
  real<lower=0> sigma_scale; // scale parameter for IG prior on sigma2
  real mu_alpha;      // prior mean for intercept (mean of y.train)
  real<lower=0> lambda_alpha;  // std dev for intercept prior (from BART)
  real<lower=0> lambda_beta;   // std dev for coefficient priors (from BART)
}

parameters {
  real alpha;               // intercept
  vector[P] beta;           // regression coefficients
  real<lower=0> sigma2;     // variance parameter
}

model {
  vector[N] mu;
  real sigma;

  // Priors (BART-style)
  alpha ~ normal(mu_alpha, lambda_alpha);  // alpha ~ N(mu_alpha, lambda_alpha^2)
  beta ~ normal(0, lambda_beta);    // beta_j ~ N(0, lambda_beta^2)

  sigma2 ~ inv_gamma(3.0/2.0, 3.0*sigma_scale/2.0); // BART-style prior

  // Derived parameter
  sigma = sqrt(sigma2);

  // Linear predictor
  mu = alpha + X * beta;

  // Likelihood for AFT model
  // power-prior weighted: RWD rows enter with weight wt
  for (n in 1:N) {
    if (event[n] == 1) {
      // Observed event
      target += wt[n] * normal_lpdf(y[n] | mu[n], sigma);
    } else {
      // Censored observation
      target += wt[n] * normal_lccdf(y[n] | mu[n], sigma);
    }
  }
}
"

# Compile Stan model once (shared across all rwd_w values)
aft_model <- stan_model(model_code = aft_stan_code)

set.seed(6)
seed <- sample(1:10000, niter*4, replace = F)

for (rwd_w in rwd_w_vals) {
  cat(sprintf("\n========== Running AFTv2 with rwd_w = %s ==========\n", rwd_w))
  ee <- 1   # reset so every rwd_w uses the same per-iteration seeds

  # Auto-load calibrated threshold when running under alternative
  threshold <- 0.95
  if (hypo == "alternative") {
    if (sc == 1 | sc == 2) {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, "_w",rwd_w,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_w",rwd_w,"_threshold.RData")
      }
    } else {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, "_cor", cor, "_w",rwd_w,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_cor", cor, "_w",rwd_w,"_threshold.RData")
      }
    }
    if (file.exists(threshold_file)) threshold <- readRDS(threshold_file)
  }

res <- data.frame(
  dist =  "lognormal",
  c = cor
)
res <- res[rep(1:nrow(res), each = niter), ]
res$bias <- NA
res$sd <- NA
res$rmst_tau <- NA
res$rmst_true <- NA
res$rmst_hat <- NA
res$bias_rmst <- NA
res$sd_rmst <- NA
res$ci_rmst <- NA
res$coverage_rmst <- NA
res$decision_rmst <- NA
res$tp_rmst <- NA
res$rmse_rmst <- NA
res$w1distance_rmst <- NA
res$w2distance_rmst <- NA
res$tp_calibrated_rmst <- NA
res$bias.trt.rmst.pop <- NA
res$sd.trt.rmst.pop <- NA
res$w2distance.trt.rmst.pop <- NA
res$bias.ctrl.rmst.pop <- NA
res$sd.ctrl.rmst.pop <- NA
res$w2distance.ctrl.rmst.pop <- NA
res$bias_subj_rmst <- NA
res$pehe_subj_rmst <- NA
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
res$bias.trt.median.pop <- NA
res$sd.trt.median.pop <- NA
res$w2distance.trt.median.pop <- NA
res$bias.ctrl.median.pop <- NA
res$sd.ctrl.median.pop <- NA
res$w2distance.ctrl.median.pop <- NA
res$iteration <- rep(1:niter, times = nrow(res) / niter)
res$rwd_w <- rwd_w   # power-prior weight on the RWD likelihood (for this run)

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
      filename_rct <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_",hypo,"_",iter,".RData")
    } else {
      filename_rct <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
    }
  } else {
    if (is.null(subsc)) {
      filename_rct <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
    } else {
      filename_rct <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
    }
  }

  # Construct filename for RWD data (same index as RCT)
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_",hypo,"_",iter,".RData")
    } else {
      filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
    }
  } else {
    if (is.null(subsc)) {
      filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
    } else {
      filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-realcov/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
    }
  }

  # Read RCT data from current iteration
  data_rct_full <- readRDS(filename_rct)
  rct_idx <- data_rct_full$X[,"D"] == 1 & data_rct_full$X[,"Z"] == 1  # single-arm: keep only RCT treated

  # Read RWD data from iteration 1
  data_rwd_full <- readRDS(filename_rwd)
  rwd_idx <- data_rwd_full$X[,"D"] == 0

  # Combine RCT (from current iter) and RWD (from iter 1)
  data$X <- rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
                  data_rwd_full$X[rwd_idx, paste0("X",1:10)])
  data$D <- c(data_rct_full$X[rct_idx, "D"],
              data_rwd_full$X[rwd_idx, "D"])
  data$Z <- c(data_rct_full$X[rct_idx, "Z"],
              data_rwd_full$X[rwd_idx, "Z"])
  data$Y <- c(data_rct_full$y[rct_idx],
              data_rwd_full$y[rwd_idx])
  data$event <- c(data_rct_full$delta[rct_idx],
                  data_rwd_full$delta[rwd_idx])

  # Center X variables based on ALL data (RCT trt + RCT ctrl + RWD ctrl)
  x_means <- colMeans(as.matrix(data$X[data$D == 1, ]))

  #------------------------------------------------
  #----------- RCT + Historical Control -----------
  #------------------------------------------------
  x_ctrl <- as.matrix(data$X[data$Z == 0, ])
  x_ctrl <- scale(x_ctrl, center = x_means, scale = FALSE)
  y_ctrl <- log(data$Y[data$Z == 0])
  event_ctrl <- data$event[data$Z == 0]

  # Estimate sigma using AFT model (BART approach)
  nu <- 3
  sigquant <- 0.90
  qchi <- qchisq(1 - sigquant, nu)

  aft_fit_ctrl <- survreg(Surv(exp(y_ctrl), event_ctrl) ~ .,
                          data = data.frame(x_ctrl),
                          dist = "lognormal")
  sigest_ctrl <- aft_fit_ctrl$scale
  lambda_ctrl <- (sigest_ctrl^2 * qchi) / nu

  # Calculate prior variances (consistent with BART)
  k <- 2.0
  total_var_ctrl <- (max(y_ctrl) - min(y_ctrl))^2 / (2 * k)^2

  # Distribute variance among parameters
  # Option 1: Equal weights (each parameter gets equal share)
  # var_per_param_ctrl <- total_var_ctrl / (ncol(x_ctrl) + 1)
  # lambda_alpha_ctrl <- sqrt(var_per_param_ctrl)
  # lambda_beta_ctrl <- sqrt(var_per_param_ctrl)

  # Option 2: More weight on beta (alpha weight = 0.5, each beta weight = 1.0)
  # w_alpha <- 0.5
  # w_beta <- 1.0
  # total_weight_ctrl <- w_alpha + ncol(x_ctrl) * w_beta
  # var_alpha_ctrl <- total_var_ctrl * w_alpha / total_weight_ctrl
  # var_beta_ctrl <- total_var_ctrl * w_beta / total_weight_ctrl
  # lambda_alpha_ctrl <- sqrt(var_alpha_ctrl)
  # lambda_beta_ctrl <- sqrt(var_beta_ctrl)

  # Option 3: More weight on alpha (alpha weight = 2.0, each beta weight = 1.0)
  w_alpha <- 2.0
  w_beta <- 1.0
  total_weight_ctrl <- w_alpha + ncol(x_ctrl) * w_beta
  var_alpha_ctrl <- total_var_ctrl * w_alpha / total_weight_ctrl
  var_beta_ctrl <- total_var_ctrl * w_beta / total_weight_ctrl
  lambda_alpha_ctrl <- sqrt(var_alpha_ctrl)
  lambda_beta_ctrl <- sqrt(var_beta_ctrl)

  # Prepare data for Stan
  stan_data_ctrl <- list(
    N = length(y_ctrl),
    P = ncol(x_ctrl),
    X = x_ctrl,
    y = y_ctrl,
    event = event_ctrl,
    wt = ifelse(data$D[data$Z == 0] == 0, rwd_w, 1.0),  # downweight RWD rows (D==0)
    sigma_scale = lambda_ctrl,  # Pass scale parameter for IG prior
    mu_alpha = 0,  # Prior mean for alpha
    lambda_alpha = lambda_alpha_ctrl,  # BART-style prior std dev for alpha
    lambda_beta = lambda_beta_ctrl     # BART-style prior std dev for beta
  )

  # Fit AFT model for control group
  fit_ctrl <- sampling(aft_model,
                       data = stan_data_ctrl,
                       chains = n_chains,
                       iter = n_iter,
                       warmup = n_warmup,
                       cores = n_cores,
                       seed = this_seed,
                       refresh = 0)

  # Extract posterior samples
  samples_ctrl <- rstan::extract(fit_ctrl)

  # Predict for target population (D=1)
  x_test <- as.matrix(data$X[data$D == 1, ])
  x_test <- scale(x_test, center = x_means, scale = FALSE)

  n_samples <- length(samples_ctrl$alpha)
  n_test <- nrow(x_test)

  # Generate predictions
  mu_pred_ctrl <- matrix(NA, n_samples, n_test)
  for (i in 1:n_samples) {
    mu_pred_ctrl[i, ] <- samples_ctrl$alpha[i] + x_test %*% samples_ctrl$beta[i, ]
  }

  # S_mix_fun and `lower` are defined here and reused by the treatment-arm root
  # finding below. The CONTROL population median is computed AFTER the treatment
  # fit, using the TREATMENT model's residual SD -- see the pop_median_draws()
  # call just before the ATE (so both arms share the RCT-scale sigma).
  S_mix_fun <- function(t, mu_vec, sig, w) {
    pnorm(log(t), mean = mu_vec, sd = sig, lower.tail = FALSE) |> as.numeric() |>
      {\(v) sum(w * v)}()
  }
  lower <- 0.05   # finite floor for the median root-search (shared with trt arm)

  #------------------------------------------------
  #---------------- RCT Treatment -----------------
  #------------------------------------------------
  x_trt <- as.matrix(data$X[data$Z == 1, ])
  x_trt <- scale(x_trt, center = x_means, scale = FALSE)
  y_trt <- log(data$Y[data$Z == 1])
  event_trt <- data$event[data$Z == 1]

  # Estimate sigma using AFT model (BART approach)
  aft_fit_trt <- survreg(Surv(exp(y_trt), event_trt) ~ .,
                         data = data.frame(x_trt),
                         dist = "lognormal")
  sigest_trt <- aft_fit_trt$scale
  lambda_trt <- (sigest_trt^2 * qchi) / nu

  # Calculate prior variances (consistent with BART)
  # Total variance for the mean function in BART
  total_var_trt <- (max(y_trt) - min(y_trt))^2 / (2 * k)^2

  # Distribute variance among parameters
  # Option 1: Equal weights (each parameter gets equal share)
  # var_per_param_trt <- total_var_trt / (ncol(x_trt) + 1)
  # lambda_alpha_trt <- sqrt(var_per_param_trt)
  # lambda_beta_trt <- sqrt(var_per_param_trt)

  # Option 2: More weight on beta (alpha weight = 0.5, each beta weight = 1.0)
  # w_alpha <- 0.5
  # w_beta <- 1.0
  # total_weight_trt <- w_alpha + ncol(x_trt) * w_beta
  # var_alpha_trt <- total_var_trt * w_alpha / total_weight_trt
  # var_beta_trt <- total_var_trt * w_beta / total_weight_trt
  # lambda_alpha_trt <- sqrt(var_alpha_trt)
  # lambda_beta_trt <- sqrt(var_beta_trt)

  # Option 3: More weight on alpha (alpha weight = 2.0, each beta weight = 1.0)
  total_weight_trt <- w_alpha + ncol(x_trt) * w_beta
  var_alpha_trt <- total_var_trt * w_alpha / total_weight_trt
  var_beta_trt <- total_var_trt * w_beta / total_weight_trt
  lambda_alpha_trt <- sqrt(var_alpha_trt)
  lambda_beta_trt <- sqrt(var_beta_trt)

  # Prepare data for Stan
  stan_data_trt <- list(
    N = length(y_trt),
    P = ncol(x_trt),
    X = x_trt,
    y = y_trt,
    event = event_trt,
    wt = rep(1.0, length(y_trt)),  # RCT treatment arm: full weight
    sigma_scale = lambda_trt,  # Pass scale parameter for IG prior
    mu_alpha = 0,  # Prior mean for alpha
    lambda_alpha = lambda_alpha_trt,  # BART-style prior std dev for alpha
    lambda_beta = lambda_beta_trt     # BART-style prior std dev for beta
  )

  # Fit AFT model for treatment group
  fit_trt <- sampling(aft_model,
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

  # Calculate median survival times for treatment using mixture root finding
  sig_draw_trt <- sqrt(samples_trt$sigma2)
  w <- rep(1/ncol(mu_pred_trt), ncol(mu_pred_trt))

  pop_median_root_trt <- rep(NA_real_, nrow(mu_pred_trt))
  upper_default <- pmin(exp(apply(mu_pred_trt, 1, max) + 6 * sig_draw_trt), 30)

  for (m in 1:nrow(mu_pred_trt)) {
    f <- function(t) S_mix_fun(t, mu_pred_trt[m, ], sig_draw_trt[m], w) - 0.5
    if (f(upper_default[m]) > 0) {
      warning(paste0("Treatment group: f(upper_default[", m, "]) > 0, clamping median to upper bound (cap)"))
      pop_median_root_trt[m] <- upper_default[m]   # clamp degenerate draw to cap (was NA)
      next
    }
    lo <- max(1e-10, lower)
    fm_lo <- f(lo); fm_hi <- f(upper_default[m])
    if (fm_lo < 0) {
      warning(paste0("Treatment group: f(lo) < 0 for m=", m, ", setting median to lower bound"))
      pop_median_root_trt[m] <- lo
    } else {
      pop_median_root_trt[m] <- uniroot(f, interval = c(lo, upper_default[m]))$root
    }
  }

  median_survival_trt <- pop_median_root_trt

  #-------------------------------------
  #---------- Calculate ATE ------------
  #-------------------------------------
  # Control counterfactual uses the TREATMENT model's residual SD (RCT-scale
  # variance), not the borrowed RWD variance, so both arms share the same sigma.
  median_survival_ctrl <- pop_median_draws(mu_pred_ctrl, sig_draw_trt)
  post_samples <- median_survival_trt / median_survival_ctrl    # n_draws
  data_tmp <- readRDS(filename_rct)
  eff <- exp(data_tmp$treat_eff_true)
  eff_star <- exp(data_tmp$treat_eff_star)

  # ---- RMST(tau)-ratio (additional output; median-ratio left unchanged) ----
  .rm <- compute_rmst_metrics(mu_pred_trt, sig_draw_trt,
                              mu_pred_ctrl, sqrt(samples_ctrl$sigma2),
                              data_tmp, tau = rmst_tau, threshold = threshold,
                              control_sigma = rmst_control_sigma)
  .sel <- which(res$c==c & res$iteration==iter)
  for (.nm in names(.rm)) res[.sel, .nm] <- .rm[[.nm]]

  # Posterior estimate
  delta_hat <- median(post_samples, na.rm = TRUE)

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
      "coverage"] <- coverage
  
  res[which(res$c==c & res$iteration==iter),
      "ci"] <- ci
  
  res[which(res$c==c & res$iteration==iter),
      "tp_calibrated"] <- tp_calibrated
  
  res[which(res$c==c & res$iteration==iter),
      "tp"] <- tp
  
  res[which(res$c==c & res$iteration==iter),
      "fp"] <- fp
  
  #------------------------------------------
  #------ Calculate subject-wise ------------
  #------------------------------------------
  post_samples_i <- exp(mu_pred_trt) / exp(mu_pred_ctrl)   # n_draws x n_subjects
  eff_i <- exp(readRDS(filename_rct)$eff_i)

  # Posterior estimates per subject
  delta_hat_i <- apply(post_samples_i, 2, median)    # n_subjects
  
  # post_samples_1i <- exp(mu_pred_trt)      # n_draws x n_subjects
  # post_samples_0i <- exp(mu_pred_ctrl)     # n_draws x n_subjects 
  # eff_i <- exp(readRDS(filename_rct)$eff_i)
  # 
  # # Posterior estimates per subject
  # delta_hat_1i <- apply(post_samples_1i, 2, median)    # n_subjects
  # delta_hat_0i <- apply(post_samples_0i, 2, median)    # n_subjects
  # delta_hat_i <- delta_hat_1i / delta_hat_0i           # n_subjects

  # Within replicate metrics
  bias_subj <- mean(delta_hat_i - eff_i)
  pehe_subj <- sqrt(mean((delta_hat_i - eff_i)^2))

  res[which(res$c==c & res$iteration==iter), "bias_subj"] <- bias_subj
  res[which(res$c==c & res$iteration==iter), "pehe_subj"] <- pehe_subj

  #-------------------------------------
  #------ Variance estimation ----------
  #-------------------------------------
  true_sigma_trt <- true_sigma_ctrl <- readRDS(filename_rct)$sigma_rct
  
  # Treatment group
  sigma_samples_trt <- sqrt(samples_trt$sigma2)
  sigma_est_trt <- mean(sigma_samples_trt)
  sd_trt_sigma <- sd(sigma_samples_trt, na.rm = TRUE)
  w2distance_trt_sigma <- sqrt(mean((sigma_samples_trt - true_sigma_trt)^2, na.rm = TRUE))
  
  res[which(res$c==c & res$iteration==iter),
      "bias.trt.sigma"] <- sigma_est_trt - true_sigma_trt
  res[which(res$c==c & res$iteration==iter),
      "sd.trt.sigma"] <- sd_trt_sigma
  res[which(res$c==c & res$iteration==iter),
      "w2distance.trt.sigma"] <- w2distance_trt_sigma

  #--------------------------------------------------
  #------ Population median survival time -----------
  #--------------------------------------------------
  true_median_trt_pop <- readRDS(filename_rct)$true_median_trt_pop
  true_median_ctrl_pop <- readRDS(filename_rct)$true_median_ctrl_pop

  # Use already calculated population medians
  median_trt_pop_hat <- median(median_survival_trt, na.rm = TRUE)
  median_ctrl_pop_hat <- median(median_survival_ctrl, na.rm = TRUE)

  # Calculate metrics
  bias_trt_median_pop <- median_trt_pop_hat - true_median_trt_pop
  sd_trt_median_pop <- sqrt(var(median_survival_trt, na.rm = TRUE))   # posterior SD (matches ATE sd)
  w2distance_trt_median_pop <- sqrt(mean((median_survival_trt - true_median_trt_pop)^2, na.rm = TRUE))
  bias_ctrl_median_pop <- median_ctrl_pop_hat - true_median_ctrl_pop
  sd_ctrl_median_pop <- sqrt(var(median_survival_ctrl, na.rm = TRUE))   # posterior SD (matches ATE sd)
  w2distance_ctrl_median_pop <- sqrt(mean((median_survival_ctrl - true_median_ctrl_pop)^2, na.rm = TRUE))

  res[which(res$c==c & res$iteration==iter),
      "bias.trt.median.pop"] <- bias_trt_median_pop
  res[which(res$c==c & res$iteration==iter),
      "sd.trt.median.pop"] <- sd_trt_median_pop
  res[which(res$c==c & res$iteration==iter),
      "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
  res[which(res$c==c & res$iteration==iter),
      "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
  res[which(res$c==c & res$iteration==iter),
      "sd.ctrl.median.pop"] <- sd_ctrl_median_pop
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

    candidate_thresholds <- seq(0.1, 0.999, by = 0.001)
    type1_errors <- sapply(candidate_thresholds, function(thresh) {
      mean(all_decisions > thresh, na.rm = TRUE)
    })

    best_idx <- which.min(abs(type1_errors - 0.05))
    calibrated_threshold_grid <- candidate_thresholds[best_idx]
    actual_type1_error <- type1_errors[best_idx]

    cat("Selected threshold: ", calibrated_threshold_grid, "\n")
    cat("Actual Type I error: ", actual_type1_error, "\n")

    if (sc == 1 | sc == 2) {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, "_w",rwd_w,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_w",rwd_w,"_threshold.RData")
      }
    } else {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, "_cor", cor, "_w",rwd_w,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-survival-realcov/res/AFTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_cor", cor, "_w",rwd_w,"_threshold.RData")
      }
    }
    saveRDS(calibrated_threshold_grid, file = threshold_file)
  }


# For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv2_p",p_obs,size_suffix,"_sc",sc,"_",hypo,"_w",rwd_w,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv2_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,"_w",rwd_w,".RData"))
    }
  } else {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv2_p",p_obs,size_suffix,"_sc",sc,"_cor",cor,"_",hypo,"_w",rwd_w,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-realcov/res/AFTv2_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",cor,"_",hypo,"_w",rwd_w,".RData"))
    }
  }

} # end for (rwd_w in rwd_w_vals)
