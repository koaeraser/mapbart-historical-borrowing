rm(list = ls())
mainDir <- "/Users/oliviazhang/Desktop/"
.scpp_cache <- file.path(tempdir(), "scpp_mBART")
dir.create(.scpp_cache, showWarnings = FALSE, recursive = TRUE)
library(Rcpp)
library(RcppEigen)
library(gtools)

data_folder <- "data_v3"  # Options: "data" or "data_v3"
RCT_ctrl <- FALSE
rct_ctrl_tag <- ""

sc <- 3
subsc <- "E"
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(mainDir, "mapbart-sim-gaussian-realcov", data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(mainDir, "mapbart-sim-gaussian-realcov", data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# Strength of unmeasured confounding (defined here -- before the calibration
# load -- so the sc/cor-matched calibration file can be selected below).
if (sc == 1 | sc == 2){
  cor <- 1
}
if (sc == 3){
  cor <- 0.7
}

# Calibrated s^2_star comes from ess_local/ess_calibration.R, which saves
# res/ess_calibration_gaussian_<version>.RData with $S_grid, $ESSbar
# (B x K), $q, and $N_target.  For each N_target multiplier m, we apply
# the Pr-rule: pick the smallest s^2 in the grid with
#   Pr{ ESSbar_b <= m * N_target } >= q.
# The result is two parallel vectors: prior_vals (s^2 values used inside
# cmbart's prior) and target_Ns (the integer target N each was
# calibrated against -- used in filenames and labels in place of the
# raw s^2 value).
N_target_multipliers <- c(1.00, 0.75, 0.50)

.cal_dir   <- file.path(mainDir, "mapbart-sim-gaussian-realcov", "ess_local", "res")
.cal_files <- list.files(.cal_dir,
                           pattern = "^ess_calibration_gaussian_.*\\.RData$",
                           full.names = TRUE)
.cal_files <- .cal_files[!grepl("_checkpoint\\.RData$", .cal_files)]
## Scenario tag mirrors the suffix ess_calibration.R appends after the config tag
## (sc1E / sc2E / sc3E_cor<rho>), so the calibration matches THIS sc/cor too.
.sc_cal_tag <- if (sc == 1 | sc == 2)
  paste0("sc", sc, if (!is.null(subsc)) subsc else "") else
  paste0("sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor)

## Per-config calibrated-s^2 selection.  Each control config (ntree,k) MUST read
## the ess calibration produced under the SAME config: s^2_star depends on
## tau_theta = (max-min)/(2*k*sqrt(ntree)) and on the leaf partial-residual
## structure, so a config's s^2 is not transferable to another config.  The
## config tag mirrors ess_calibration.R's cfg_tag = sprintf("_nt%d_k%g", ntree, k)
## and is REQUIRED here (sits BEFORE the scenario suffix .sc_cal_tag).  Called
## inside the control-config loop so prior_vals/target_Ns are per-config.
pick_prior_for <- function(ctrl_ntree, ctrl_k) {
  cfg_tag   <- sprintf("nt%d_k%g", ctrl_ntree, ctrl_k)
  cal_match <- .cal_files[grepl(paste0("^ess_calibration_gaussian_.*_", sprintf("nt%d_k%g", ctrl_ntree, ctrl_k), "_", .sc_cal_tag, "\\.RData$"),
                                basename(.cal_files))]
  if (length(cal_match)) {
    cal_latest <- cal_match[order(file.info(cal_match)$mtime, decreasing = TRUE)][1]
    cal        <- readRDS(cal_latest)
    .pick <- function(m) {
      Ghat_m <- colMeans(cal$ESSbar <= m * cal$N_target, na.rm = TRUE)
      ok     <- Ghat_m >= cal$q
      if (any(ok)) cal$S_grid[which(ok)[1]] else NA_real_
    }
    pv     <- sapply(N_target_multipliers, .pick)
    tn     <- as.integer(round(N_target_multipliers * cal$N_target))
    ok_idx <- which(!is.na(pv))
    pv     <- pv[ok_idx]
    tn     <- tn[ok_idx]
    cat(sprintf("mapbart-sim-gaussian-realcov/mBART.R: [%s] from %s\n",
                cfg_tag, basename(cal_latest)))
    for (j in seq_along(pv))
      cat(sprintf("  N=%d -> prior(s^2) = %.4g\n", tn[j], pv[j]))
  } else {
    # Fallback when no config-/sc-matched ess calibration exists.
    pv <- 0.05
    tn <- 100L
    cat(sprintf(paste0("mapbart-sim-gaussian-realcov/mBART.R: [%s] WARNING -- no config-/sc-matched ",
                "calibration (ess_calibration_gaussian_*_%s_%s.RData) in %s.\n  Using fallback ",
                "s^2 = %.4g (target_N=%d) + the 1e-4 (N=Inf) reference.  Run ess_calibration ",
                "with --ntree %d --k %g first for calibrated targets.\n"),
                cfg_tag, cfg_tag, .sc_cal_tag, .cal_dir, pv, tn, ctrl_ntree, ctrl_k))
  }
  # Also run the near-zero discrepancy reference (s^2 = 0.0001, i.e. maximal
  # borrowing -> N = Inf) alongside the calibrated s^2_star values.  The
  # reference run is tagged _NInf_; the calibrated runs keep their integer-N tags.
  pv <- c(pv, 0.0001)
  tn <- c(tn, Inf)
  cat(sprintf("mapbart-sim-gaussian-realcov/mBART.R: [%s] s^2 sweep (calibrated s^2_star + 0.0001) = %s\n",
              cfg_tag, paste(pv, collapse=", ")))
  list(prior_vals = pv, target_Ns = tn)
}

threshold_default <- 0.95

# Infer niter from the number of simulated data files matching this scenario
data_prefix <- if (sc == 1 | sc == 2) {
  paste0("data_p", p_obs, "_sc", sc, if (!is.null(subsc)) subsc else "", "_", hypo, "_")
} else {
  paste0("data_p", p_obs, "_sc", sc, if (!is.null(subsc)) subsc else "", "_cor", cor, "_", hypo, "_")
}
niter <- length(Sys.glob(file.path(mainDir, "mapbart-sim-gaussian-realcov", data_folder, paste0(data_prefix, "*.RData"))))
if (niter == 0) stop("No simulated data files found for this scenario")

ndpost=1500L
nskip=100L
keepevery=1L

# alpha_beta <- data.frame(alpha = c(0.9, 0.95, 0.98),  # -> bigger trees
#                          beta = c(2.5, 2, 1.5))       # -> bigger trees
# ntree <- c(1, 5, 10, 50, 100, 200)
alpha_beta <- data.frame(alpha = c(0.95),  # -> bigger trees
                         beta = c(2))       # -> bigger trees
ntree <- c(50)   # treatment-arm BART ntree (control arm uses the per-config ctrl_ntree below)

# Control-arm BART config sweep (ntree:k pairs), mirroring mapbart-case-study-mm/mBART_analysis.R.
# Each config flexes/shrinks the CONTROL-arm BART differently (smaller ntree + larger k
# => more shrinkage) and writes its own result files tagged _nt<ntree>_k<k>.  The
# treatment arm is unaffected.  run_all.R overrides `mbart_control_configs`.
mbart_control_configs <- "5:2,10:2,50:2"   # comma-separated ntree:k pairs (control-arm BART)
control_configs <- lapply(strsplit(mbart_control_configs, "[, ]+")[[1]], function(p) {
  if (!nzchar(p)) return(NULL)
  v <- strsplit(p, ":")[[1]]
  if (length(v) != 2L || anyNA(suppressWarnings(as.numeric(v))))
    stop("Bad mbart_control_configs entry '", p, "' (expected ntree:k, e.g. 50:2)")
  list(ntree = as.integer(v[1]), k = as.numeric(v[2]))
})
control_configs <- Filter(Negate(is.null), control_configs)
stopifnot(length(control_configs) >= 1)

# Only run configs that have a config-/sc-matched ess calibration file (same
# required tag as pick_prior_for); skip the rest rather than running them on the
# uncalibrated fallback s^2.
has_cal_for <- function(ctrl_ntree, ctrl_k) {
  cfg_tag <- sprintf("nt%d_k%g", ctrl_ntree, ctrl_k)
  any(grepl(paste0("^ess_calibration_gaussian_.*_", cfg_tag, "_", .sc_cal_tag, "\\.RData$"),
            basename(.cal_files)))
}
.has_cal <- vapply(control_configs, function(cf) has_cal_for(cf$ntree, cf$k), logical(1))
if (any(!.has_cal))
  cat(sprintf("mBART.R: skipping control configs with no %s calibration: %s\n",
              .sc_cal_tag,
              paste(sapply(control_configs[!.has_cal],
                           function(cf) sprintf("(ntree=%d,k=%g)", cf$ntree, cf$k)), collapse = ", ")))
control_configs <- control_configs[.has_cal]
if (length(control_configs) < 1) {
  # No config calibrated for THIS scenario -> skip MAP-BART for this sc entirely
  # (the control-config loop below is a no-op), rather than hard-stopping a
  # multi-scenario run_all sweep.
  warning(sprintf("mBART.R: no control config has a matching %s ess calibration (ess_calibration_gaussian_*_nt<ntree>_k<k>_%s.RData in %s) -- SKIPPING MAP-BART for this scenario. Run ess_calibration for this sc/config first.",
                  .sc_cal_tag, .sc_cal_tag, .cal_dir))
} else {
  cat(sprintf("mBART.R: control-arm configs = %s\n",
              paste(sapply(control_configs, function(cf) sprintf("(ntree=%d,k=%g)", cf$ntree, cf$k)),
                    collapse = ", ")))
}

# Single-arm (RCT_ctrl=FALSE): use mBART.real, the external-control-only build of
# cmbart (same signature), mirroring mapbart-case-study-mm/mBART_analysis.R. Otherwise use mBART.
if (isTRUE(RCT_ctrl)) {
  sourceCpp(paste0(mainDir,"/mBART/cmbart.cpp"), cacheDir = .scpp_cache)
} else {
  sourceCpp(paste0(mainDir,"/mBART.real/cmbart.cpp"), cacheDir = .scpp_cache)
}
sourceCpp(paste0(mainDir,"/wBART/cwbart.cpp"), cacheDir = .scpp_cache)

for (.cfg in control_configs) {       # NEW: sweep control-arm BART configs
ctrl_ntree <- .cfg$ntree
ctrl_k     <- .cfg$k
cat(sprintf("\n######## CONTROL-ARM CONFIG: ntree=%d, k=%g ########\n", ctrl_ntree, ctrl_k))

# Read THIS config's calibrated s^2 sweep (config-matched; see pick_prior_for).
.pp        <- pick_prior_for(ctrl_ntree, ctrl_k)
prior_vals <- .pp$prior_vals
target_Ns  <- .pp$target_Ns

# Auto-load calibrated threshold when running under alternative (config-tagged;
# the threshold file does not depend on prior in the Gaussian real variant).
threshold <- threshold_default
if (hypo == "alternative") {
  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/mBART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, rct_ctrl_tag,"_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/mBART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, subsc, rct_ctrl_tag,"_threshold.RData")
    }
  } else {
    if (is.null(subsc)) {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/mBART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, "_cor", cor, rct_ctrl_tag,"_threshold.RData")
    } else {
      threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/mBART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, subsc, "_cor", cor, rct_ctrl_tag,"_threshold.RData")
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
              filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
            } else {
              filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
            }
          } else {
            if (is.null(subsc)) {
              filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
            } else {
              filename_rct <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
            }
          }
          
          # Construct filename for RWD data (current iteration)
          if (sc == 1 | sc == 2) {
            if (is.null(subsc)) {
              filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_",hypo,"_",iter,".RData")
            } else {
              filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_",hypo,"_",iter,".RData")
            }
          } else {
            if (is.null(subsc)) {
              filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter,".RData")
            } else {
              filename_rwd <- paste0(mainDir,"/mapbart-sim-gaussian-realcov/",data_folder,"/data_p",p_obs,"_sc",sc,subsc,"_cor",c,"_",hypo,"_",iter,".RData")
            }
          }
          
          # Read RCT data from current iteration
          data_rct_full <- readRDS(filename_rct)
          rct_idx <- data_rct_full$X[,"D"] == 1
          if (!RCT_ctrl) rct_idx <- rct_idx & data_rct_full$X[,"Z"] == 1  # single-arm: keep only RCT treated
          
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
          event <- as.integer(rep(1, length(data$Y)))
          
          # Also need filename variable for reading true values later
          filename <- filename_rct
          
          #------------------------------------------------
          #----------- RCT + Historical Control -----------
          #------------------------------------------------
          type = 1
          x.train = data.frame(X = as.matrix(data$X[data$Z == 0, ]))
          s.train = data$D[data$Z == 0]
          y.train = data$Y[data$Z == 0]
          event_train = as.integer(event[data$Z == 0])
          
          x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))
          s.test = rep(1, nrow(x.test))
          
          xinfo=matrix(0,0,0)
          usequants=FALSE
          rm.const=TRUE
          
          sigdf=3
          sigquant=.90
          power=2
          base=0.95
          
          n = length(y.train)
          
          numcut=100L
          transposed=FALSE
          
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
          
          offset0 = 0
          offset1 = 1
          y.train = (y.train - offset0) / offset1
          
          # RWD group (s==0)
          df0 = data.frame(t(x.train[, s.train == 0]))
          df0$y.train = y.train[s.train == 0]
          lmf0 = lm(y.train ~ ., df0)
          sigest0 = summary(lmf0)$sigma

          # RCT group (s==1). Single-arm (RCT_ctrl=FALSE) has no concurrent
          # control, so the s==1 group is empty; cmbart (mBART.real) ignores the
          # s==1 prior in that case, so fall back to the RWD estimate to keep the
          # passed sigest1/sigscale1 finite.
          if (sum(s.train == 1) > 0) {
            df1 = data.frame(t(x.train[, s.train == 1]))
            df1$y.train = y.train[s.train == 1]
            lmf1 = lm(y.train ~ ., df1)
            sigest1 = summary(lmf1)$sigma
          } else {
            sigest1 = sigest0
          }
          
          # Calculate separate sigscale for each group
          sigscale1 = (sigest1^2) * qchisq(0.10, df = sigdf) / sigdf # P(\sigma^2 < sigest^2) = 0.9
          sigscale0 = (sigest0^2) * qchisq(0.10, df = sigdf) / sigdf # P(\sigma^2 < sigest^2) = 0.9
          
          t=ctrl_k
          tau = (max(y.train)-min(y.train))/(2*t*sqrt(ctrl_ntree))

          set.seed(this_seed)
          res_ctrl = cmbart(type,
                            n,
                            p,
                            np,
                            x.train,
                            y.train,
                            s.train,
                            event_train,
                            x.test,
                            s.test,
                            ctrl_ntree,
                            numcut,
                            ndpost*keepevery,
                            nskip,
                            keepevery,
                            power,
                            base,
                            offset0,
                            offset1,
                            sigdf,
                            sigscale1,
                            sigscale0,
                            sigest1,
                            sigest0,
                            tau,
                            xinfo,
                            "IG",  # sigma_prior: "IG" for inverse gamma, "hC" for half-Cauchy
                            prior,  # prior parameter for rate(3.0*prior/2)
                            seed[ee]  # seed for random number generation
          )

          # mBART.real (single-arm) returns only the shared residual SD as
          # $sigma_rwd (the concurrent-control variance is unidentified with no
          # internal control). Alias it to $sigma_rct so the downstream thinning
          # and sigma metrics work -- matches mapbart-case-study-mm/mBART_analysis.R, which uses
          # the shared sigma_rwd in the no-concurrent-control case.
          if (!isTRUE(RCT_ctrl)) res_ctrl$sigma_rct <- res_ctrl$sigma_rwd

          if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
          if(keepevery>1) res_ctrl$sigma_rct = c(res_ctrl$sigma_rct[nskip.],
                                                 res_ctrl$sigma_rct[nskip+seq(1, ndpost*keepevery, keepevery)])
          res_ctrl$sigma_rct = res_ctrl$sigma_rct[-(nskip.)] 
          
          #------------------------------------------------
          #---------------- RCT Treatment -----------------
          #------------------------------------------------
          x.train = data.frame(X = as.matrix(data$X[data$Z == 1, ]))
          y.train = data$Y[data$Z == 1]
          x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))
          
          sparse=FALSE
          theta=0
          omega=1
          a=0.5
          b=1
          augment=FALSE
          rho=NULL
          xinfo=matrix(0.0,0,0)
          usequants=FALSE
          cont=FALSE
          rm.const=TRUE
          sigest=NA
          sigdf=3
          sigquant=.90
          # k=2.0
          power=beta
          base=alpha
          k=2.0
          sigmaf=NA
          lambda=NA
          fmean=0
          w=rep(1,length(y.train))
          numcut=100L
          
          nkeeptrain=ndpost
          nkeeptest=ndpost
          nkeeptestmean=ndpost
          nkeeptreedraws=ndpost
          printevery=100L
          transposed=FALSE
          
          n = length(y.train)
          
          if(!transposed) {
            source(paste0(mainDir,"/bartModelMatrix.R"))
            temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                                   cont=cont, xinfo=xinfo, rm.const=rm.const)
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
          
          y.train = y.train-fmean
          if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
            nkeeptrain=ndpost
          }
          if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
            nkeeptest=ndpost
          }
          if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
            nkeeptestmean=ndpost
          }
          if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
            nkeeptreedraws=ndpost
          }
          #--------------------------------------------------
          nu=sigdf
          if(is.na(lambda)) {
            if(is.na(sigest)) {
              if(p < n) {
                df = data.frame(t(x.train),y.train)
                lmf = lm(y.train~.,df)
                sigest = summary(lmf)$sigma
              } else {
                sigest = sd(y.train)
              }
            }
            qchi = qchisq(1.0-sigquant,nu)
            lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
          } else {
            sigest=sqrt(lambda)
          }
          
          if(is.na(sigmaf)) {
            tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
          } else {
            tau = sigmaf/sqrt(ntree)
          }
          
          
          set.seed(this_seed)
          res_trt = cwbart(n,  #number of observations in training data
                           p,  #dimension of x
                           np, #number of observations in test data
                           x.train,   #pxn training data x
                           y.train,   #pxn training data x
                           x.test,   #p*np test data x
                           ntree,
                           numcut,
                           ndpost*keepevery,
                           nskip,
                           power,
                           base,
                           tau,
                           nu,
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
                           nkeeptrain,
                           nkeeptest,
                           nkeeptestmean,
                           nkeeptreedraws,
                           printevery,
                           xinfo
          )
          
          if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
          if(keepevery>1) res_trt$sigma = c(res_trt$sigma[nskip.],
                                            res_trt$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
          res_trt$sigma = res_trt$sigma[-(nskip.)] 
          
          #-------------------------------------
          #---------- Calculate ATE ------------
          #-------------------------------------
          post_samples <- rowMeans(res_trt$yhat.test - res_ctrl$f.test)
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
              "ci"] <- ci
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "coverage"] <- coverage
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "tp_calibrated"] <- tp_calibrated
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "tp"] <- tp
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "fp"] <- fp
          
          #------------------------------------------
          #------ Calculate subject-wise ------------
          #------------------------------------------
          post_samples_i <- res_trt$yhat.test - res_ctrl$f.test   # n_draws x n_subjects
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
          
          # Treatment group
          sigma_samples_trt <- res_trt$sigma
          sigma_est_trt <- mean(sigma_samples_trt)
          sd_trt_sigma <- sd(sigma_samples_trt)
          w2distance_trt_sigma <- sqrt(mean((sigma_samples_trt - true_sigma_trt)^2, na.rm = TRUE))
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "bias.trt.sigma"] <- sigma_est_trt - true_sigma_trt
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "sd.trt.sigma"] <- sd_trt_sigma
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "w2distance.trt.sigma"] <- w2distance_trt_sigma
          
          # Control group
          sigma_samples_ctrl <- res_ctrl$sigma_rct
          sigma_est_ctrl <- mean(sigma_samples_ctrl)
          sd_ctrl_sigma <- sd(sigma_samples_ctrl)
          w2distance_ctrl_sigma <- sqrt(mean((sigma_samples_ctrl - true_sigma_trt)^2, na.rm = TRUE))
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "bias.ctrl.sigma"] <- sigma_est_ctrl - true_sigma_ctrl
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "sd.ctrl.sigma"] <- sd_ctrl_sigma
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "w2distance.ctrl.sigma"] <- w2distance_ctrl_sigma
          
          #--------------------------------------------------
          #------ Population mean estimation ----------------
          #--------------------------------------------------
          true_mean_trt_pop <- readRDS(filename)$true_mean_trt_pop
          true_mean_ctrl_pop <- readRDS(filename)$true_mean_ctrl_pop
          
          # Calculate population means
          mean_trt_pop_hat <- mean(rowMeans(res_trt$yhat.test))
          mean_ctrl_pop_hat <- mean(rowMeans(res_ctrl$f.test))
          
          # Calculate metrics
          mu_pred_trt_samples <- rowMeans(res_trt$yhat.test)
          mu_pred_ctrl_samples <- rowMeans(res_ctrl$f.test)
          bias_trt_median_pop <- mean_trt_pop_hat - true_mean_trt_pop
          w2distance_trt_median_pop <- sqrt(mean((mu_pred_trt_samples - true_mean_trt_pop)^2))
          bias_ctrl_median_pop <- mean_ctrl_pop_hat - true_mean_ctrl_pop
          w2distance_ctrl_median_pop <- sqrt(mean((mu_pred_ctrl_samples - true_mean_ctrl_pop)^2))
          
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "bias.trt.median.pop"] <- bias_trt_median_pop
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
          res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==iter),
              "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop
          
          print(paste0("Done for iteration ", iter, " cor = ", c))
        }
        
        print(paste0("Done for cor = ", c,
                     " H = ", H,
                     " alpha = ", alpha,
                     " beta = ", beta))
      }}}
  
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
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, rct_ctrl_tag,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, subsc, rct_ctrl_tag,"_threshold.RData")
      }
    } else {
      if (is.null(subsc)) {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, "_cor", c, rct_ctrl_tag,"_threshold.RData")
      } else {
        threshold_file <- paste0(mainDir, "/mapbart-sim-gaussian-realcov/res/MAP-BART_p", p_obs, "_nt", ctrl_ntree, "_k", ctrl_k, "_sc", sc, subsc, "_cor", c, rct_ctrl_tag,"_threshold.RData")
      }
    }
    saveRDS(calibrated_threshold_grid, file = threshold_file)
  }
  
  # For sc == 1 or 2: no cor suffix; for sc == 3: include cor suffix
  if (sc == 1 | sc == 2) {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,"_nt",ctrl_ntree,"_k",ctrl_k,"_sc",sc,"_N",target_N,"_",hypo,rct_ctrl_tag,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,"_nt",ctrl_ntree,"_k",ctrl_k,"_sc",sc,subsc,"_N",target_N,"_",hypo,rct_ctrl_tag,".RData"))
    }
  } else {
    if (is.null(subsc)) {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,"_nt",ctrl_ntree,"_k",ctrl_k,"_sc",sc,"_cor",c,"_N",target_N,"_",hypo,rct_ctrl_tag,".RData"))
    } else {
      saveRDS(res, file=paste0(mainDir,"/mapbart-sim-gaussian-realcov/res/MAP-BART_p",p_obs,"_nt",ctrl_ntree,"_k",ctrl_k,"_sc",sc,subsc,"_cor",c,"_N",target_N,"_",hypo,rct_ctrl_tag,".RData"))
    }
  }

} # End of (prior, target_N) loop

} # NEW: End of control-arm config loop
