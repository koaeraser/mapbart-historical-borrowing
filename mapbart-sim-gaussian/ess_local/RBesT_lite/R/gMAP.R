## Trimmed gMAP for the single-call signature used in simulate_ess.R:
##
##   gMAP(formula = cbind(y.bar, y.se) ~ 1 | study,
##        data, beta.prior = <scalar sd>,
##        tau.dist = "InvGamma", tau.prior = c(a, b),
##        family = gaussian,
##        version = c("standard", "hybrid", "synthetic"))
##
## Three model versions are supported (Stan source in ../stan/):
##
##   "standard"  -- the default MAP prior:
##                  theta_h = mu + tau * xi_h, mu ~ N(0, beta.prior^2)
##                  theta_star ~ N(mu, tau^2)
##
##   "hybrid"    -- two-arm hybrid design (mu fixed at 0):
##                  theta_h    = tau * xi_h
##                  theta_star ~ N(0, tau^2)
##                  beta.prior is ignored.
##
##   "synthetic" -- single-arm synthetic control (single anchor):
##                  theta_h    ~ N(0, beta.prior^2)
##                  theta_star ~ N(theta_h, tau^2)
##                  no exchangeable random-effect machinery.
##
## tau ~ InvGamma(a, b) in all three cases.  Other RBesT options
## (binomial/poisson, other tau dists, weights, prior_PD, etc.) are
## removed.  The returned list contains only the fields downstream
## code (mixfit.gMAP) reads: $family, $thin, $fit, $sigma_ref.
##
## REPRODUCIBILITY: rstan::sampling() is given an explicit `seed=`
## drawn from R's RNG state, so once R's seed has been set upstream
## each gMAP() call is fully deterministic.

gMAP <- function(formula, data,
                 beta.prior, tau.prior,
                 tau.dist = "InvGamma",
                 family   = gaussian,
                 version  = c("standard", "hybrid", "synthetic"),
                 iter     = getOption("RBesT.MC.iter",   6000),
                 warmup   = getOption("RBesT.MC.warmup", 2000),
                 thin     = getOption("RBesT.MC.thin",   4),
                 init     = getOption("RBesT.MC.init",   1),
                 chains   = getOption("RBesT.MC.chains", 4),
                 cores    = getOption("mc.cores",        1L)) {

  version <- match.arg(version)
  stopifnot(identical(tau.dist, "InvGamma"))
  fam <- if (is.function(family)) family() else family
  stopifnot(fam$family == "gaussian", fam$link == "identity")

  ## --- model frame ----------------------------------------------------
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  f  <- Formula::Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())

  ## response is a 2-col matrix (y.bar, y.se)
  y_mat <- model.response(mf)
  y     <- array(y_mat[, 1])
  y_se  <- array(y_mat[, 2])
  H     <- NROW(y)

  ## grouping factor (RHS after the bar)
  group.factor <- factor(model.part(f, data = mf, rhs = 2)[, 1])
  group.index  <- array(as.integer(group.factor))
  n.groups     <- nlevels(group.factor)

  ## --- design matrix (intercept) -------------------------------------
  X  <- model.matrix(f, mf, rhs = 1)
  mX <- NCOL(X)

  ## --- guesses used by Stan rescaling --------------------------------
  ## sigma_ref / sigma_guess: pooled sd of the supplied SEs
  sigma_ref   <- sqrt(H / sum(1 / y_se^2))
  sigma_guess <- sigma_ref
  tau_guess   <- if (n.groups > 1)
                   max(sigma_guess / 10, sd(tapply(y, group.index, mean)))
                 else 1

  ## --- beta prior (scalar sd, mean 0) --------------------------------
  if (!missing(beta.prior) && length(beta.prior) == 1)
    beta.prior <- matrix(c(0, beta.prior), nrow = 1,
                         dimnames = list(NULL, c("mean", "sd")))
  if (!missing(beta.prior) && !is.matrix(beta.prior))
    beta.prior <- matrix(beta.prior, mX, 2, byrow = TRUE,
                         dimnames = list(NULL, c("mean", "sd")))

  ## --- tau prior (InvGamma) -----------------------------------------
  tau.prior <- as.numeric(tau.prior)
  stopifnot(length(tau.prior) == 2)

  ## --- pooled fit (used to seed the Stan rescale init for default) --
  fit.pooled <- glm.fit(X, y, weights = as.vector(1 / y_se^2),
                        offset = rep(0, H), family = fam)

  ## --- rescaling guesses for tau and beta ----------------------------
  nInf <- 0.9 * (sigma_guess / tau_guess)^2
  if (n.groups > 1) {
    ms <- square_root_gamma_stats((n.groups - 1) / 2,
                                  2 * tau_guess^2 / (n.groups - 1))
    ms[2] <- sqrt(1 + 2 / (n.groups - 1)) * ms[2]
    tau_raw_guess <- c(log(ms[1]) - log(sqrt(1 + ms[2]^2 / ms[1]^2)),
                       sqrt(log(1 + ms[2]^2 / ms[1]^2)))
  } else {
    tau_raw_guess <- c(log(tau_guess), 1)
  }
  beta_raw_guess <- rbind(mean = fit.pooled$coefficients,
                          sd   = rep(sigma_guess / sqrt(nInf), mX))

  rescale <- getOption("RBesT.MC.rescale", TRUE)
  if (!rescale) {
    tau_raw_guess[2]   <- 1
    beta_raw_guess[2,] <- 1
  }

  ## --- assemble Stan data -- one shape per version ------------------
  dataL <- switch(version,
    standard  = list(H = H, mX = mX, X = X,
                   y = y, y_se = y_se,
                   group_index = group.index, n_groups = n.groups,
                   beta_prior = beta.prior, tau_prior = tau.prior,
                   tau_raw_guess  = tau_raw_guess,
                   beta_raw_guess = beta_raw_guess),
    hybrid    = list(H = H,
                   y = y, y_se = y_se,
                   group_index = group.index, n_groups = n.groups,
                   tau_prior = tau.prior,
                   tau_raw_guess = tau_raw_guess),
    synthetic = list(H = H,
                   y = y, y_se = y_se,
                   beta_prior = as.numeric(beta.prior[1, ]),
                   tau_prior = tau.prior,
                   tau_raw_guess = tau_raw_guess)
  )

  control <- modifyList(list(adapt_delta = 0.99, stepsize = 0.01,
                             max_treedepth = 20),
                        getOption("RBesT.MC.control", list()))

  ## --- Stan ---------------------------------------------------------
  ## Pull a Stan seed from R's current RNG state so a single
  ## set.seed() call upstream determines both R-side randomness AND
  ## Stan's NUTS trajectory.  Without this, rstan::sampling() draws
  ## its own seed from a non-RNG source on some platforms.
  stan_seed <- sample.int(.Machine$integer.max, 1L)

  exclude_pars <- switch(version,
                         standard  = c("beta_raw", "tau_raw", "xi_eta"),
                         hybrid    = c("tau_raw", "xi_eta"),
                         synthetic = c("tau_raw"))
  capture.output(
    fit <- rstan::sampling(
      stanmodels[[version]], data = dataL,
      seed = stan_seed,
      warmup = warmup, iter = iter, chains = chains, cores = cores,
      thin = thin, init = init, refresh = 0, control = control,
      algorithm = "NUTS", open_progress = FALSE,
      pars = exclude_pars, include = FALSE,
      save_warmup = getOption("RBesT.MC.save_warmup", FALSE)
    )
  )
  if (attributes(fit)$mode != 0) stop("Stan sampler did not run successfully!")

  ## mixfit.gMAP only reads $family, $thin, $fit, $sigma_ref.
  structure(list(family = fam, thin = 1, fit = fit,
                 sigma_ref = sigma_ref, version = version),
            class = "gMAP")
}

## Helper used by the rescale guesses above.
square_root_gamma_stats <- function(a, b) {
  m <- sqrt(b) * exp(lgamma(0.5 + a) - lgamma(a))
  v <- b * a - m^2
  c(mean = m, sd = sqrt(v))
}
