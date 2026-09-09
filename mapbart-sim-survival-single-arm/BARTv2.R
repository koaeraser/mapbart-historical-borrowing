rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop/"
.scpp_cache <- file.path(tempdir(), "scpp_BARTv2")
dir.create(.scpp_cache, showWarnings = FALSE, recursive = TRUE)
library(Rcpp)
library(RcppEigen)
library(gtools)
library(survival)
library(dplyr)

source(paste0(mainDir, "/mapbart-sim-survival-single-arm/rmst_helpers.R"))
rmst_tau <- 5   # RMST restriction horizon (admin-censoring horizon); additional output
rmst_control_sigma <- "own"  # control-arm RMST sigma: "trt"|"own"|"adaptive" (own = control's own sigma-hat; avoids inflated trt sigma at small n)

data_folder <- "data_v3"  # Options: "data" or "data_v3"
size_option <- "default"   # "default" (RCT trt 200) or "small" (RCT trt 30)
size_suffix <- if (size_option == "small") "_n30" else "_n200"

sc <- 1
subsc <- "E"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-survival-single-arm", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-survival-single-arm", data_folder))
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

# Infer niter from the number of simulated data files matching this scenario
data_prefix <- if (sc == 1 | sc == 2) {
  paste0("data_p", p_obs, size_suffix, "_sc", sc, if (!is.null(subsc)) subsc else "", "_", hypo, "_")
} else {
  paste0("data_p", p_obs, size_suffix, "_sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor, "_", hypo, "_")
}
niter <- length(Sys.glob(file.path(mainDir, "mapbart-sim-survival-single-arm", data_folder, paste0(data_prefix, "*.RData"))))
if (niter == 0) stop("No simulated data files found for this scenario")

# Auto-load calibrated threshold when running under alternative
if (hypo == "alternative") {
  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, "_cor", cor, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_cor", cor, "_threshold.RData")
    }
  }

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}

ndpost=1500L
nskip=100L
keepevery=1L

# alpha_beta <- data.frame(alpha = c(0.9, 0.95, 0.98),  # -> bigger trees
#                          beta = c(2.5, 2, 1.5))       # -> bigger trees
# ntree <- c(1, 5, 10, 50, 100, 200)
alpha_beta <- data.frame(alpha = c(0.95),  # -> bigger trees
                         beta = c(2))       # -> bigger trees
ntree <- c(50)

sourceCpp(paste0(mainDir,"/aBART/cabart.cpp"), cacheDir = .scpp_cache)

set.seed(6)
seed <- sample(1:10000, niter*4, replace = F)
ee <- 1

res <- expand.grid(
  id = 1:nrow(alpha_beta),
  H = ntree,
  c = cor
)

# Merge with alpha-beta combinations
res <- merge(res,
             alpha_beta,
             by.x = "id",
             by.y = "row.names",
             all.x = TRUE)
res <- res[ , c("alpha", "beta", "H", "c")]
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

# For threshold calibration under H0
if (hypo == "null") {
  all_decisions <- numeric(nrow(alpha_beta) * length(cor) * length(ntree) * niter)
  decision_idx <- 1
}

cat("Loaded threshold from H0 run:", threshold)

for (which in 1:nrow(alpha_beta)){
  alpha <- alpha_beta[which, "alpha"]
  beta <- alpha_beta[which, "beta"]
  for (c in cor){
    for (H in ntree){
      for (iter in 1:niter){
        this_seed <- seed[ee]   # one seed per iteration; both arms share it
        set.seed(this_seed)
        ee <- ee + 1
        
        data = list()
        
        # Construct filename for RCT data (current iteration)
        # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
        if (sc == 1 | sc == 2) {
          if (is.null(subsc)) {
            filename_rct <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_",hypo,"_",iter,".RData")
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
          }
        } else {
          if (is.null(subsc)) {
            filename_rct <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
          } else {
            filename_rct <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
          }
        }

        # Construct filename for RWD data (same index as RCT)
        if (sc == 1 | sc == 2) {
          if (is.null(subsc)) {
            filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_",hypo,"_",iter,".RData")
          } else {
            filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
          }
        } else {
          if (is.null(subsc)) {
            filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
          } else {
            filename_rwd <- paste0(mainDir,"/mapbart-sim-survival-single-arm/",data_folder,"/data_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
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
        
        #------------------------------------------------
        #----------- RCT + Historical Control -----------
        #------------------------------------------------
        x.train = data.frame(X = as.matrix(data$X[data$Z == 0, ]))
        y.train = log(data$Y[data$Z == 0])
        event = as.integer(data$event[data$Z == 0])

        x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))
        
        ntype=1
        sparse=FALSE
        theta=0
        omega=1
        a=0.5
        b=1
        augment=FALSE
        rho=NULL
        xinfo=matrix(0.0,0,0)
        usequants=FALSE
        rm.const=TRUE
        sigdf=3
        sigquant=.90
        k=2.0
        power=beta
        base=alpha
        sigmaf=NA
        lambda=NA
        offset=0
        w=rep(1,length(y.train))
        numcut=100L
        printevery=10000L
        transposed=FALSE
        
        n = length(y.train)
        
        if(!transposed) {
          source(paste0(mainDir,"/bartModelMatrix.R"))
          temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                                 xinfo=xinfo, rm.const=rm.const)
          x.train = t(temp$X)
          numcut = temp$numcut
          xinfo = temp$xinfo
          if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
          }
          rm.const <- temp$rm.const
          grp <- temp$grp
          rm(temp)
        }
        
        p = nrow(x.train)
        np = ncol(x.test)
        if(length(rho)==0) rho=p
        if(length(rm.const)==0) rm.const <- 1:p
        if(length(grp)==0) grp <- 1:p
        
        y.train = y.train-offset
        if(is.na(lambda)) {
          # sigest = summary(lm(y.train~.,
          #                     data.frame(t(x.train),y.train)))$sigma
          aft_fit = survreg(Surv(exp(y.train), event) ~ .,
                            data = data.frame(t(x.train)),
                            dist = "lognormal")
          sigest = aft_fit$scale
          
          qchi = qchisq(1.0-sigquant,sigdf)
          lambda = (sigest*sigest*qchi)/sigdf #lambda parameter for sigma prior
        } else {
          sigest=sqrt(lambda)
        }
        
        if(is.na(sigmaf)) {
          tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
          tau = sigmaf/sqrt(ntree)
        }
        
        set.seed(this_seed)
        res_ctrl = cabart(ntype,
                          n,
                          p,
                          np,
                          x.train,
                          y.train,
                          event,
                          x.test,
                          ntree,
                          numcut,
                          ndpost*keepevery,
                          nskip,
                          keepevery,
                          power,
                          base,
                          offset,
                          tau,
                          sigdf,
                          lambda,
                          sigest,
                          w,
                          sparse,
                          theta,
                          omega,
                          grp,
                          a,
                          b,
                          rho,
                          augment,
                          printevery,
                          xinfo,
                          seed[ee]  # seed for random number generation
        )

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
        x.train = data.frame(X = as.matrix(data$X[data$Z == 1, ]))
        y.train = log(data$Y[data$Z == 1])
        event = as.integer(data$event[data$Z == 1])
        x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))
        
        ntype=1
        sparse=FALSE
        theta=0
        omega=1
        a=0.5
        b=1
        augment=FALSE
        rho=NULL
        xinfo=matrix(0.0,0,0)
        usequants=FALSE
        rm.const=TRUE
        sigest=NA
        sigdf=3
        sigquant=.90
        k=2.0
        power=beta
        base=alpha
        sigmaf=NA
        lambda=NA
        offset=0
        w=rep(1,length(y.train))
        numcut=100L
        printevery=10000L
        transposed=FALSE
        
        n = length(y.train)
        
        if(!transposed) {
          source(paste0(mainDir,"/bartModelMatrix.R"))
          temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                                 xinfo=xinfo, rm.const=rm.const)
          x.train = t(temp$X)
          numcut = temp$numcut
          xinfo = temp$xinfo
          if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
          }
          rm.const <- temp$rm.const
          grp <- temp$grp
          rm(temp)
        }
        
        p = nrow(x.train)
        np = ncol(x.test)
        if(length(rho)==0) rho=p
        if(length(rm.const)==0) rm.const <- 1:p
        if(length(grp)==0) grp <- 1:p
        
        y.train = y.train-offset
        
        if(is.na(lambda)) {
          # sigest = summary(lm(y.train~.,
          #                     data.frame(t(x.train),y.train)))$sigma
          aft_fit = survreg(Surv(exp(y.train), event) ~ .,
                            data = data.frame(t(x.train)),
                            dist = "lognormal")
          sigest = aft_fit$scale
          
          qchi = qchisq(1.0-sigquant,sigdf)
          lambda = (sigest*sigest*qchi)/sigdf #lambda parameter for sigma prior
        } else {
          sigest=sqrt(lambda)
        }
        
        if(is.na(sigmaf)) {
          tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
          tau = sigmaf/sqrt(ntree)
        }
        
        
        set.seed(this_seed)
        res_trt = cabart(ntype,
                         n,
                         p,
                         np,
                         x.train,
                         y.train,
                         event,
                         x.test,
                         ntree,
                         numcut,
                         ndpost*keepevery,
                         nskip,
                         keepevery,
                         power,
                         base,
                         offset,
                         tau,
                         sigdf,
                         lambda,
                         sigest,
                         w,
                         sparse,
                         theta,
                         omega,
                         grp,
                         a,
                         b,
                         rho,
                         augment,
                         printevery,
                         xinfo,
                         seed[ee]  # seed for random number generation
        )
        
        if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
        if(keepevery>1) res_trt$sigma = c(res_trt$sigma[nskip.],
                                          res_trt$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
        res_trt$sigma = res_trt$sigma[-(nskip.)]
        
        # Calculate median survival times for treatment
        mu_draws <- res_trt$yhat.test
        sig_draw <- res_trt$sigma
        w <- rep(1/ncol(mu_draws), ncol(mu_draws))
        
        pop_median_root_trt <- rep(NA_real_, nrow(mu_draws))
        upper_default <- pmin(exp(apply(mu_draws, 1, max) + 6 * sig_draw), 30)
        
        for (m in 1:nrow(mu_draws)) {
          f <- function(t) S_mix_fun(t, mu_draws[m, ], sig_draw[m], w) - 0.5
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
        median_survival_ctrl <- pop_median_draws(res_ctrl$yhat.test, res_trt$sigma)
        post_samples <- median_survival_trt / median_survival_ctrl
        data_tmp <- readRDS(filename_rct)
        eff <- exp(data_tmp$treat_eff_true)
        eff_star <- exp(data_tmp$treat_eff_star)

        # ---- RMST(tau)-ratio (additional output; median-ratio left unchanged) ----
        .rm <- compute_rmst_metrics(res_trt$yhat.test, res_trt$sigma,
                                    res_ctrl$yhat.test, tail(res_ctrl$sigma, length(res_trt$sigma)),
                                    data_tmp, tau = rmst_tau, threshold = threshold,
                                    control_sigma = rmst_control_sigma)
        .sel <- which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter)
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
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "bias"] <- bias
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "sd"] <- SD
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "rmse"] <- rmse
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "w1distance"] <- w1distance
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "w2distance"] <- w2distance
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "coverage"] <- coverage
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "ci"] <- ci
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "tp_calibrated"] <- tp_calibrated
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "tp"] <- tp
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "fp"] <- fp
        
        #------------------------------------------
        #------ Calculate subject-wise ------------
        #------------------------------------------
        post_samples_i <- exp(res_trt$yhat.test) / exp(res_ctrl$yhat.test)   # n_draws x n_subjects
        eff_i <- exp(readRDS(filename_rct)$eff_i)

        # Posterior estimates per subject
        delta_hat_i <- apply(post_samples_i, 2, median)    # n_subjects
        
        # post_samples_1i <- exp(res_trt$yhat.test)      # n_draws x n_subjects
        # post_samples_0i <- exp(res_ctrl$yhat.test)     # n_draws x n_subjects 
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
        sigma_samples_trt <- res_trt$sigma
        sigma_est_trt <- mean(sigma_samples_trt)
        sd_trt_sigma <- sd(sigma_samples_trt, na.rm = TRUE)
        w2distance_trt_sigma <- sqrt(mean((sigma_samples_trt - true_sigma_trt)^2, na.rm = TRUE))
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "bias.trt.sigma"] <- sigma_est_trt - true_sigma_trt
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "sd.trt.sigma"] <- sd_trt_sigma
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
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

        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "bias.trt.median.pop"] <- bias_trt_median_pop
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "sd.trt.median.pop"] <- sd_trt_median_pop
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "sd.ctrl.median.pop"] <- sd_ctrl_median_pop
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
            "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop

        print(paste0("Done for iteration ", iter, " cor = ", c))
      }
      print(paste0("Done for cor = ", c,
                   " H = ", H,
                   " alpha = ", alpha,
                   " beta = ", beta))
    }
  }
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
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, "_cor", c, "_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-survival-single-arm/res/BARTv2_p", p_obs, size_suffix, "_sc", sc, subsc, "_cor", c, "_threshold.RData")
    }
  }

  saveRDS(calibrated_threshold_grid, file = threshold_file)
}


# For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
if (sc == 1 | sc == 2) {
  if (is.null(subsc)) {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-single-arm/res/BARTv2_p",p_obs,size_suffix,"_sc",sc,"_",hypo,".RData"))
  } else {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-single-arm/res/BARTv2_p",p_obs,size_suffix,"_sc",sc,subsc,"_",hypo,".RData"))
  }
} else {
  if (is.null(subsc)) {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-single-arm/res/BARTv2_p",p_obs,size_suffix,"_sc",sc,"_cor",c,"_",hypo,".RData"))
  } else {
    saveRDS(res, file=paste0(mainDir,"/mapbart-sim-survival-single-arm/res/BARTv2_p",p_obs,size_suffix,"_sc",sc,subsc,"_cor",c,"_",hypo,".RData"))
  }
}

