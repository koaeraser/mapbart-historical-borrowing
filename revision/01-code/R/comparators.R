# comparators.R: the paper's comparator methods with a common interface.
#
# Interface. Every method is fit_<name>(dat, outcome, ...) where dat is the
# list returned by gen_data() and outcome is "gaussian" or "survival". It
# returns the object built by make_fit() (estimands.R):
#   $draws    data.frame with B rows: ate, trt, ctrl (Gaussian) or ate_med,
#             trt_med, ctrl_med, ate_rmst, trt_rmst, ctrl_rmst (survival);
#   $surface  posterior mean of the control conditional mean at all n_rct RCT
#             profiles (treated first, then control), outcome scale;
#   $map      optional P(borrow) at each RCT control profile (NULL here; the
#             proposed method fills it);
#   $ess      optional list(prior = , realized = ).
# Methods are registered in METHODS by name; register_method() adds one, so
# the proposed LRC-BART can be plugged in without touching the harness.
#
# Common settings: B = 1000 post-warm-up draws after 1000 burn-in for every
# sampler (the paper gives no explicit count). All Gaussian linear models and
# the AFT models are fitted by exact conjugate Gibbs samplers (the LM/AFT
# posteriors are conjugate given sigma^2, and censoring is handled by data
# augmentation, so the sampler targets the exact posterior). The paper says
# "all Bayesian methods were fitted by Hamiltonian Monte Carlo"; the target
# posterior is identical, only the algorithm differs. stan/hier_lm.stan holds
# the HMC version of HierLM for cross-checking.
#
# Deviations from the paper, and points it leaves unspecified (all listed):
#  1. HierLM/HierAFT: the paper writes beta_{g,j} ~ N(mu_beta, tau_j^2) with a
#     single mu_beta; we use one hypermean per coefficient, mu_{beta,j} ~
#     N(0, 10^2), which is the natural reading of a coefficient-wise
#     hierarchy. The residual variances sigma_g^2 receive the BART-style
#     Inv-Gamma(nu/2, nu*lambda_g/2) prior with lambda_g from a group-specific
#     OLS (or survreg) fit, as the paper states for HierAFT.
#  2. AFT-NP/CP/PP priors: the paper gives alpha ~ N(0, tau_alpha^2), beta ~
#     N(0, tau_beta^2 I) without values; we use the LM values (5^2 and 2^2).
#  3. MAP prior: the paper states tau ~ Half-Normal(0, 1) and three targets
#     N_target in {100, 75, 50}. With one external study of n = 300 a fixed
#     HN(0, 1) scale gives an ELIR prior ESS of about 2, which cannot produce
#     the borrowing seen in the paper's MAP rows. We therefore calibrate the
#     half-normal scale so that the ELIR ESS of the mixture-approximated MAP
#     prior equals N_target (falling back to the largest attainable ESS when
#     the target is unreachable). The mixture, ESS, and posterior update use
#     RBesT (mixnorm, ess, postmix). The treatment arm and the RWD summary are
#     unadjusted sample means, as the paper describes.
#  4. PSCL is implemented from Wang et al. (2020) as read through the paper:
#     propensity score for D by logistic regression on X1..X10 over RCT and
#     RWD subjects, trimming of RWD outside the RCT propensity range, 5
#     strata by RCT propensity quintiles, the total borrowed size A = 100
#     split across strata in proportion to the stratum overlap coefficient of
#     the propensity distributions (the "distance" allocation), power
#     parameter lambda_s = min(A_s, m_s) / m_s, composite-likelihood stratum
#     means and their sandwich variance, and standardization by the RCT
#     stratum shares. It is a frequentist estimator; to fit the interface its
#     "draws" are N(estimate, SE^2) samples, and the arm-specific draws are
#     provided although the paper reports none for PSCL.
#  5. BART uses BART::wbart (Gaussian) and BART::abart (log-normal AFT) with
#     H = 50 trees, alpha = 0.95, beta = 2, k = 2, nu = 3, q = 0.90 as in the
#     paper. The treatment arm always gets its own ensemble on the RCT
#     treated subjects; the pooled ensembles (CP, PP) use the pooled residual
#     scale for the control arm.

suppressPackageStartupMessages({
  library(BART)
  library(survival)
  library(RBesT)
})

N_BURN <- 1000L
N_DRAW <- 1000L
XVARS <- paste0("X", 1:10)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

xmat <- function(df, vars = XVARS) as.matrix(df[, vars, drop = FALSE])

# Inverse-gamma draw for sigma^2 ~ IG(a, b).
rinvgamma1 <- function(a, b) 1 / rgamma(1, shape = a, rate = b)

# BART-style calibration of the residual-variance prior scale: lambda such
# that P(sigma < sigma_hat) = q under sigma^2 ~ IG(nu/2, nu*lambda/2), i.e.
# lambda = sigma_hat^2 * qchisq(1 - q, nu) / nu.
lambda_bart <- function(sigma_hat, nu = 3, q = 0.90) sigma_hat^2 * qchisq(1 - q, nu) / nu

# Preliminary residual scale: OLS for Gaussian, survreg(lognormal) for survival.
sigma_prelim <- function(y, X, status = NULL) {
  if (is.null(status)) {
    fit <- lm.fit(cbind(1, X), y)
    return(sqrt(sum(fit$residuals^2) / max(1, length(y) - ncol(X) - 1)))
  }
  fit <- try(survreg(Surv(exp(y), status) ~ X, dist = "lognormal"), silent = TRUE)
  if (inherits(fit, "try-error")) return(sd(y))
  fit$scale
}

# Draw log survival times for censored observations from N(mu, sigma^2)
# truncated below at the observed log censoring time a.
rtnorm_lower <- function(mu, sigma, a) {
  p_up <- pnorm((a - mu) / sigma, lower.tail = FALSE)
  u <- runif(length(mu)) * p_up
  mu + sigma * qnorm(u, lower.tail = FALSE)
}

# Design and response for one arm/source set. For survival, y = log(time).
prep_xy <- function(df, outcome, add_d = FALSE) {
  X <- xmat(df)
  if (add_d) X <- cbind(X, D = df$d)
  if (outcome == "gaussian") {
    list(X = X, y = df$y, status = NULL)
  } else {
    list(X = X, y = log(df$time), status = df$status)
  }
}

# ---------------------------------------------------------------------------
# Bayesian linear regression / log-normal AFT by conjugate Gibbs
# ---------------------------------------------------------------------------
# y = alpha + X beta + eps, eps ~ N(0, sigma^2). Priors alpha ~ N(0, sd_a^2),
# beta ~ N(0, sd_b^2 I), sigma^2 ~ IG(nu/2, nu*lambda/2) with lambda from the
# preliminary fit. Censored y (status = 0) are augmented. Returns draws of
# theta = (alpha, beta) (B x (p+1)) and sigma (B).
gibbs_lm <- function(y, X, status = NULL, sd_a = 5, sd_b = 2, nu = 3, q = 0.90,
                     n_burn = N_BURN, n_draw = N_DRAW) {
  n <- length(y); p <- ncol(X)
  Xd <- cbind(1, X)
  prior_prec <- c(1 / sd_a^2, rep(1 / sd_b^2, p))
  sig_hat <- sigma_prelim(y, X, status)
  lambda <- lambda_bart(sig_hat, nu, q)
  cens <- if (is.null(status)) rep(FALSE, n) else status == 0
  y_aug <- y
  XtX <- crossprod(Xd)
  # initial values
  theta <- qr.solve(Xd, y)
  sigma2 <- sig_hat^2
  out_theta <- matrix(NA_real_, n_draw, p + 1)
  out_sigma <- numeric(n_draw)
  for (it in seq_len(n_burn + n_draw)) {
    if (any(cens)) {
      mu_c <- as.vector(Xd[cens, , drop = FALSE] %*% theta)
      y_aug[cens] <- rtnorm_lower(mu_c, sqrt(sigma2), y[cens])
    }
    V_inv <- XtX / sigma2 + diag(prior_prec)
    R <- chol(V_inv)
    m <- backsolve(R, backsolve(R, crossprod(Xd, y_aug) / sigma2, transpose = TRUE))
    theta <- as.vector(m + backsolve(R, rnorm(p + 1)))
    ssr <- sum((y_aug - Xd %*% theta)^2)
    sigma2 <- rinvgamma1((nu + n) / 2, (nu * lambda + ssr) / 2)
    if (it > n_burn) {
      out_theta[it - n_burn, ] <- theta
      out_sigma[it - n_burn] <- sqrt(sigma2)
    }
  }
  list(theta = out_theta, sigma = out_sigma)
}

# Predicted conditional means at test profiles (B x n_test).
predict_lm <- function(fit, X_test) fit$theta %*% t(cbind(1, X_test))

# The treatment-arm model shared by the linear-family methods: fitted on the
# RCT treated subjects, evaluated at all RCT profiles.
fit_trt_lm <- function(dat, outcome) {
  d <- prep_xy(dat$rct_trt, outcome)
  f <- gibbs_lm(d$y, d$X, d$status)
  list(mu = predict_lm(f, dat$x_rct), sigma = f$sigma)
}

fit_linear_family <- function(dat, outcome, pooling = c("NP", "CP", "PP")) {
  pooling <- match.arg(pooling)
  trt <- fit_trt_lm(dat, outcome)
  ctrl_df <- switch(pooling,
                    NP = dat$rct_ctrl,
                    CP = rbind(dat$rct_ctrl, dat$rwd_ctrl),
                    PP = rbind(dat$rct_ctrl, dat$rwd_ctrl))
  d <- prep_xy(ctrl_df, outcome, add_d = (pooling == "PP"))
  f <- gibbs_lm(d$y, d$X, d$status)
  X_test <- if (pooling == "PP") cbind(dat$x_rct, D = 1) else dat$x_rct
  mu0 <- predict_lm(f, X_test)
  make_fit(trt$mu, mu0, outcome, trt$sigma, f$sigma, dat$constants$tau)
}

fit_lm_np <- function(dat, outcome, ...) fit_linear_family(dat, outcome, "NP")
fit_lm_cp <- function(dat, outcome, ...) fit_linear_family(dat, outcome, "CP")
fit_lm_pp <- function(dat, outcome, ...) fit_linear_family(dat, outcome, "PP")

# ---------------------------------------------------------------------------
# Hierarchical linear / AFT model (HierLM, HierAFT) by conjugate Gibbs
# ---------------------------------------------------------------------------
# y_i = alpha_g + X_i beta_g + eps_i, eps_i ~ N(0, sigma_g^2), g = 1 (RCT
# control), 2 (RWD). theta_g = (alpha_g, beta_g) ~ N(mu, diag(tau^2)),
# mu_k ~ N(0, 10^2), tau_k^2 ~ IG(3/2, 3*0.5/2), sigma_g^2 ~ IG(nu/2,
# nu*lambda_g/2). Returns draws of theta_1 (RCT control) and sigma_1.
gibbs_hier <- function(y1, X1, y2, X2, status1 = NULL, status2 = NULL,
                       nu = 3, q = 0.90, hyper_sd = 10, tau_a = 1.5, tau_b = 0.75,
                       n_burn = N_BURN, n_draw = N_DRAW) {
  p <- ncol(X1); K <- p + 1
  ys <- list(y1, y2); Xs <- list(cbind(1, X1), cbind(1, X2))
  st <- list(status1, status2)
  ns <- c(length(y1), length(y2))
  cens <- lapply(1:2, function(g) if (is.null(st[[g]])) rep(FALSE, ns[g]) else st[[g]] == 0)
  y_aug <- ys
  XtX <- lapply(Xs, crossprod)
  sig_hat <- c(sigma_prelim(y1, X1, status1), sigma_prelim(y2, X2, status2))
  lambda <- lambda_bart(sig_hat, nu, q)
  theta <- lapply(1:2, function(g) qr.solve(Xs[[g]], ys[[g]]))
  sigma2 <- sig_hat^2
  mu <- (theta[[1]] + theta[[2]]) / 2
  tau2 <- rep(0.5, K)
  out_theta <- matrix(NA_real_, n_draw, K)
  out_sigma <- numeric(n_draw)
  for (it in seq_len(n_burn + n_draw)) {
    for (g in 1:2) {
      if (any(cens[[g]])) {
        mu_c <- as.vector(Xs[[g]][cens[[g]], , drop = FALSE] %*% theta[[g]])
        y_aug[[g]][cens[[g]]] <- rtnorm_lower(mu_c, sqrt(sigma2[g]), ys[[g]][cens[[g]]])
      }
      V_inv <- XtX[[g]] / sigma2[g] + diag(1 / tau2)
      R <- chol(V_inv)
      rhs <- crossprod(Xs[[g]], y_aug[[g]]) / sigma2[g] + mu / tau2
      m <- backsolve(R, backsolve(R, rhs, transpose = TRUE))
      theta[[g]] <- as.vector(m + backsolve(R, rnorm(K)))
      ssr <- sum((y_aug[[g]] - Xs[[g]] %*% theta[[g]])^2)
      sigma2[g] <- rinvgamma1((nu + ns[g]) / 2, (nu * lambda[g] + ssr) / 2)
    }
    # hypermeans and hypervariances
    prec <- 2 / tau2 + 1 / hyper_sd^2
    mean_mu <- ((theta[[1]] + theta[[2]]) / tau2) / prec
    mu <- rnorm(K, mean_mu, sqrt(1 / prec))
    ss <- (theta[[1]] - mu)^2 + (theta[[2]] - mu)^2
    tau2 <- 1 / rgamma(K, shape = tau_a + 1, rate = tau_b + ss / 2)
    if (it > n_burn) {
      out_theta[it - n_burn, ] <- theta[[1]]
      out_sigma[it - n_burn] <- sqrt(sigma2[1])
    }
  }
  list(theta = out_theta, sigma = out_sigma)
}

fit_hier <- function(dat, outcome, ...) {
  trt <- fit_trt_lm(dat, outcome)
  d1 <- prep_xy(dat$rct_ctrl, outcome)
  d2 <- prep_xy(dat$rwd_ctrl, outcome)
  f <- gibbs_hier(d1$y, d1$X, d2$y, d2$X, d1$status, d2$status)
  mu0 <- predict_lm(f, dat$x_rct)
  make_fit(trt$mu, mu0, outcome, trt$sigma, f$sigma, dat$constants$tau)
}

# ---------------------------------------------------------------------------
# BART family (wbart for Gaussian, abart for log-normal AFT)
# ---------------------------------------------------------------------------
BART_OPTS <- list(ntree = 50L, k = 2, power = 2, base = 0.95, sigdf = 3,
                  sigquant = 0.90, ndpost = N_DRAW, nskip = N_BURN)

bart_fit_one <- function(df, outcome, X_test, add_d = FALSE) {
  X <- xmat(df)
  if (add_d) X <- cbind(X, D = df$d)
  # Suppress the package's progress printing.
  if (outcome == "gaussian") {
    fit <- suppressMessages(capture_bart(
      wbart(x.train = X, y.train = df$y, x.test = X_test,
            ntree = BART_OPTS$ntree, k = BART_OPTS$k, power = BART_OPTS$power,
            base = BART_OPTS$base, sigdf = BART_OPTS$sigdf,
            sigquant = BART_OPTS$sigquant, ndpost = BART_OPTS$ndpost,
            nskip = BART_OPTS$nskip, printevery = 100000L)))
    sig <- tail(fit$sigma, BART_OPTS$ndpost)
  } else {
    fit <- suppressMessages(capture_bart(
      abart(x.train = X, times = df$time, delta = df$status, x.test = X_test,
            K = 2, ntree = BART_OPTS$ntree, k = BART_OPTS$k, power = BART_OPTS$power,
            base = BART_OPTS$base, sigdf = BART_OPTS$sigdf,
            sigquant = BART_OPTS$sigquant, ndpost = BART_OPTS$ndpost,
            nskip = BART_OPTS$nskip, keepevery = 1L, printevery = 100000L)))
    sig <- tail(fit$sigma, BART_OPTS$ndpost)
  }
  list(mu = fit$yhat.test, sigma = sig)
}

capture_bart <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- expr))
  out
}

fit_bart_family <- function(dat, outcome, pooling = c("NP", "CP", "PP")) {
  pooling <- match.arg(pooling)
  trt <- bart_fit_one(dat$rct_trt, outcome, dat$x_rct)
  ctrl_df <- if (pooling == "NP") dat$rct_ctrl else rbind(dat$rct_ctrl, dat$rwd_ctrl)
  X_test <- if (pooling == "PP") cbind(dat$x_rct, D = 1) else dat$x_rct
  ctrl <- bart_fit_one(ctrl_df, outcome, X_test, add_d = (pooling == "PP"))
  make_fit(trt$mu, ctrl$mu, outcome, trt$sigma, ctrl$sigma, dat$constants$tau)
}

fit_bart_np <- function(dat, outcome, ...) fit_bart_family(dat, outcome, "NP")
fit_bart_cp <- function(dat, outcome, ...) fit_bart_family(dat, outcome, "CP")
fit_bart_pp <- function(dat, outcome, ...) fit_bart_family(dat, outcome, "PP")

# ---------------------------------------------------------------------------
# MAP prior (RBesT), Gaussian only
# ---------------------------------------------------------------------------
# Mixture approximation of the MAP prior for the RCT control mean built from
# the single RWD study summary (ybar_2, se_2) and tau ~ HN(0, s): a K-point
# equal-probability quadrature over tau gives an RBesT mixnorm with
# components N(ybar_2, se_2^2 + tau_k^2). The scale s is calibrated so that
# ess(., method = "elir") equals N_target (see deviation 3 in the header).
map_prior_mix <- function(ybar2, se2, s, sigma_ref, K = 20) {
  pk <- (seq_len(K) - 0.5) / K
  tau_k <- s * qnorm(0.5 + pk / 2)
  comp <- rbind(rep(1 / K, K), rep(ybar2, K), sqrt(se2^2 + tau_k^2))
  m <- do.call(mixnorm, c(lapply(seq_len(K), function(k) comp[, k]), list(sigma = sigma_ref)))
  m
}

map_calibrate_scale <- function(ybar2, se2, sigma_ref, n_target) {
  ess_at <- function(log_s) {
    ess(map_prior_mix(ybar2, se2, exp(log_s), sigma_ref), method = "elir",
        sigma = sigma_ref) - n_target
  }
  lo <- log(1e-3); hi <- log(20)
  if (ess_at(lo) <= 0) return(exp(lo))   # target unreachable: borrow at most the study
  if (ess_at(hi) >= 0) return(exp(hi))
  exp(uniroot(ess_at, c(lo, hi), tol = 1e-3)$root)
}

fit_map <- function(dat, outcome, n_target = 100, ...) {
  if (outcome != "gaussian") stop("MAP is defined for the Gaussian study only")
  y1 <- dat$rct_ctrl$y; y2 <- dat$rwd_ctrl$y; yt <- dat$rct_trt$y
  n1 <- length(y1); n2 <- length(y2); nt <- length(yt)
  sigma_ref <- sd(y1)
  ybar2 <- mean(y2); se2 <- sd(y2) / sqrt(n2)
  s <- map_calibrate_scale(ybar2, se2, sigma_ref, n_target)
  prior <- map_prior_mix(ybar2, se2, s, sigma_ref)
  ess_prior <- ess(prior, method = "elir", sigma = sigma_ref)
  post <- postmix(prior, m = mean(y1), se = sd(y1) / sqrt(n1))
  ctrl <- rmix(post, N_DRAW)
  trt <- rnorm(N_DRAW, mean(yt), sd(yt) / sqrt(nt))
  # Marginal-mean method: the "surface" is flat at the posterior control mean.
  list(draws = data.frame(ate = trt - ctrl, trt = trt, ctrl = ctrl),
       surface = rep(mean(ctrl), nrow(dat$x_rct)),
       map = NULL, ess = list(prior = ess_prior, tau_scale = s), extra = NULL)
}

fit_map_100 <- function(dat, outcome, ...) fit_map(dat, outcome, n_target = 100)
fit_map_75 <- function(dat, outcome, ...) fit_map(dat, outcome, n_target = 75)
fit_map_50 <- function(dat, outcome, ...) fit_map(dat, outcome, n_target = 50)

# ---------------------------------------------------------------------------
# PSCL (Wang et al. 2020), Gaussian only
# ---------------------------------------------------------------------------
# Overlap coefficient of two samples by a common histogram.
overlap_coef <- function(a, b, nbins = 10) {
  if (length(a) == 0 || length(b) == 0) return(0)
  br <- seq(min(a, b), max(a, b), length.out = nbins + 1)
  br[1] <- br[1] - 1e-8; br[nbins + 1] <- br[nbins + 1] + 1e-8
  pa <- tabulate(cut(a, br, labels = FALSE), nbins) / length(a)
  pb <- tabulate(cut(b, br, labels = FALSE), nbins) / length(b)
  sum(pmin(pa, pb))
}

fit_pscl <- function(dat, outcome, A_total = 100, n_strata = 5, ...) {
  if (outcome != "gaussian") stop("PSCL is defined for the Gaussian study only")
  rct <- rbind(dat$rct_trt, dat$rct_ctrl)
  rwd <- dat$rwd_ctrl
  all <- rbind(rct, rwd)
  ps_fit <- glm(d ~ ., data = all[, c(XVARS, "d")], family = binomial())
  ps <- fitted(ps_fit)
  ps_rct <- ps[all$d == 1]; ps_rwd <- ps[all$d == 0]
  keep_rwd <- ps_rwd >= min(ps_rct) & ps_rwd <= max(ps_rct)
  rwd <- rwd[keep_rwd, ]; ps_rwd <- ps_rwd[keep_rwd]
  cuts <- quantile(ps_rct, probs = seq(0, 1, length.out = n_strata + 1))
  cuts[1] <- -Inf; cuts[n_strata + 1] <- Inf
  s_rct <- cut(ps_rct, cuts, labels = FALSE)
  s_rwd <- cut(ps_rwd, cuts, labels = FALSE)
  # Allocation of the borrowed size across strata by overlap ("distance").
  r_s <- vapply(seq_len(n_strata), function(s) overlap_coef(ps_rct[s_rct == s], ps_rwd[s_rwd == s]), 0)
  if (sum(r_s) <= 0) r_s <- rep(1, n_strata)
  A_s <- A_total * r_s / sum(r_s)
  m_s <- tabulate(s_rwd, n_strata)
  A_s <- pmin(A_s, m_s)
  lambda_s <- ifelse(m_s > 0, A_s / m_s, 0)
  # Residual scales (pooled within source).
  zr <- rct$z
  s1 <- sd(rct$y[zr == 0]); s2 <- sd(rwd$y); st <- sd(rct$y[zr == 1])
  w_s <- tabulate(s_rct, n_strata) / length(s_rct)
  th_c <- v_c <- th_t <- v_t <- numeric(n_strata)
  for (s in seq_len(n_strata)) {
    yc <- rct$y[zr == 0 & s_rct == s]; yr <- rwd$y[s_rwd == s]; yt <- rct$y[zr == 1 & s_rct == s]
    n_s <- length(yc); ms <- length(yr); l <- lambda_s[s]
    den <- n_s + l * ms
    if (den > 0) {
      th_c[s] <- (sum(yc) + l * sum(yr)) / den
      v_c[s] <- (n_s * s1^2 + l^2 * ms * s2^2) / den^2
    } else { th_c[s] <- NA; v_c[s] <- NA }
    if (length(yt) > 0) { th_t[s] <- mean(yt); v_t[s] <- st^2 / length(yt) } else { th_t[s] <- NA; v_t[s] <- NA }
  }
  ok <- !is.na(th_c) & !is.na(th_t)
  w <- w_s[ok] / sum(w_s[ok])
  est_c <- sum(w * th_c[ok]); se_c <- sqrt(sum(w^2 * v_c[ok]))
  est_t <- sum(w * th_t[ok]); se_t <- sqrt(sum(w^2 * v_t[ok]))
  ctrl <- rnorm(N_DRAW, est_c, se_c)
  trt <- rnorm(N_DRAW, est_t, se_t)
  list(draws = data.frame(ate = trt - ctrl, trt = trt, ctrl = ctrl),
       surface = rep(est_c, nrow(dat$x_rct)), map = NULL,
       ess = list(prior = sum(A_s), lambda = lambda_s), extra = NULL)
}

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
METHODS <- list(
  "LM-NP" = fit_lm_np, "LM-CP" = fit_lm_cp, "LM-PP" = fit_lm_pp,
  "BART-NP" = fit_bart_np, "BART-CP" = fit_bart_cp, "BART-PP" = fit_bart_pp,
  "PSCL" = fit_pscl,
  "MAP-100" = fit_map_100, "MAP-75" = fit_map_75, "MAP-50" = fit_map_50,
  "HierLM" = fit_hier,
  # survival aliases (same code, the outcome argument selects the likelihood)
  "AFT-NP" = fit_lm_np, "AFT-CP" = fit_lm_cp, "AFT-PP" = fit_lm_pp,
  "HierAFT" = fit_hier
)

register_method <- function(name, fun) {
  stopifnot(is.function(fun))
  METHODS[[name]] <<- fun
  invisible(name)
}

get_method <- function(name) {
  if (!name %in% names(METHODS)) stop("unknown method: ", name)
  METHODS[[name]]
}

# The paper's comparator set for each outcome type.
default_methods <- function(outcome) {
  if (outcome == "gaussian") {
    c("LM-NP", "LM-CP", "LM-PP", "BART-NP", "BART-CP", "BART-PP", "PSCL",
      "MAP-100", "MAP-75", "MAP-50", "HierLM")
  } else {
    c("AFT-NP", "AFT-CP", "AFT-PP", "BART-NP", "BART-CP", "BART-PP", "HierAFT")
  }
}
