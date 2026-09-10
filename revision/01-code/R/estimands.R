# estimands.R: turn posterior draws of the arm-specific conditional means at
# the RCT covariate profiles into draws of the paper's standardized estimands.
#
# Every method must produce, for each posterior draw b, the conditional mean
# mu_a^(b)(x_i) at every RCT profile i (treated and control profiles alike,
# n_rct of them) for both arms a = 1 (treatment) and a = 0 (control). The
# standardization to the RCT covariate distribution is then the average over
# the n_rct profiles.
#
# Gaussian:  ate = mean_i mu_1(x_i) - mean_i mu_0(x_i), trt, ctrl.
# Survival:  log-normal AFT with mu_a(x) the location and sigma_a the scale;
#            the population survival curve S_a(t) = mean_i S_a(t | x_i); the
#            estimands are the population median ratio (ate_med) and the RMST
#            ratio at horizon tau (ate_rmst), with the arm-specific medians and
#            RMSTs. The survival helpers lnorm_rmst and lnorm_pop_median live in
#            dgp.R so that the truth and the estimates use one implementation.

# mu1, mu0: B x n_rct matrices of draws. sig1, sig0: length-B vectors (survival
# only). Returns a data.frame of B rows.
standardized_draws <- function(mu1, mu0, outcome, sig1 = NULL, sig0 = NULL, tau = 3) {
  if (outcome == "gaussian") {
    trt <- rowMeans(mu1); ctrl <- rowMeans(mu0)
    return(data.frame(ate = trt - ctrl, trt = trt, ctrl = ctrl))
  }
  B <- nrow(mu1)
  stopifnot(length(sig1) == B, length(sig0) == B)
  med1 <- vapply(seq_len(B), function(b) lnorm_pop_median(mu1[b, ], sig1[b]), 0)
  med0 <- vapply(seq_len(B), function(b) lnorm_pop_median(mu0[b, ], sig0[b]), 0)
  rm1 <- rowMeans(lnorm_rmst(mu1, sig1, tau))
  rm0 <- rowMeans(lnorm_rmst(mu0, sig0, tau))
  data.frame(ate_med = med1 / med0, trt_med = med1, ctrl_med = med0,
             ate_rmst = rm1 / rm0, trt_rmst = rm1, ctrl_rmst = rm0)
}

# Assemble the common fit object. surface is the posterior mean of the control
# conditional mean at all RCT profiles (outcome scale; log scale for survival).
make_fit <- function(mu1, mu0, outcome, sig1 = NULL, sig0 = NULL, tau = 3,
                     map = NULL, ess = NULL, extra = NULL) {
  list(draws = standardized_draws(mu1, mu0, outcome, sig1, sig0, tau),
       surface = colMeans(mu0),
       map = map, ess = ess, extra = extra)
}
