# dgp.R: data-generating mechanisms for the LRC-BART revision simulation.
#
# Implements the paper's Scenarios 1 to 3 (Section "Simulation Studies" of
# paper/main.tex) exactly as written, plus the revision's Sc4 (partial
# compatibility) and Sc5 (fully shifted RWD) from 00-design/method_spec.md
# Section E. Gaussian and censored log-normal (survival) outcomes.
#
# Design choices that the paper leaves implicit (documented here so they can
# be audited):
#   * Sample sizes are exact: a pool of 2 * (n_rct + n_rwd) subjects receives
#     the source indicator D from the scenario's mechanism, and exactly n_rct
#     subjects are sampled from D = 1 and n_rwd from D = 0. This preserves the
#     covariate distribution given D and gives the paper's n_RCT = n_RWD = 300.
#   * Treatment within the RCT is Bernoulli(ratio / (ratio + 1)), so the arm
#     sizes are "approximately" 200/100 (2:1) or 225/75 (3:1) as in the paper.
#   * Under Sc3 the observed X5, X6 are the proxies rho * U_std + sqrt(1 - rho^2)
#     * eps (mean 0), and the outcome uses W = X with (X5, X6) replaced by
#     (U1, U2), as the paper states.
#   * The true arm-specific means are returned in two versions: standardized
#     to the replicate's own RCT covariate profiles ("sample", primary, this is
#     what every fitted estimand is standardized to) and to the RCT population
#     given D = 1 by Monte Carlo ("pop").
#   * Survival truths use the RCT residual scale sigma_RCT for both arms
#     (both arms are RCT arms).

DGP_MU <- c(1, -1, 0.5, -1, 2, -2, 2, -3, 1, 0)
DGP_BETA <- c(-0.50, -0.75, -0.50, -0.50, -1.35, -0.80, -0.01, -0.01, -0.01, -0.01)
DGP_BETA_D <- -1.20

# Bounded bump centred at the confounder means.
g_bump <- function(u) u^2 / (1 + (u / 1.5)^2)

# The paper's nonlinear mean function f(W; beta). W is an n x 10 matrix.
f_mean <- function(W, beta = DGP_BETA) {
  beta[1] * W[, 1]^2 +
    beta[2] * exp(-0.5 * (W[, 2] + 1)^2) +
    beta[3] * abs(W[, 3] - 0.5) +
    beta[4] * abs(W[, 4] + 1) +
    beta[5] * tanh(W[, 5] - 2) +
    beta[6] * tanh(W[, 6] + 2) +
    2 * (g_bump(W[, 5] - 2) + g_bump(W[, 6] + 2)) +
    as.vector(W[, 7:10] %*% beta[7:10])
}

# Scenario-specific constants.
dgp_constants <- function(scenario, outcome) {
  base <- if (scenario == 1 || scenario >= 4) "sc1" else "sc23"
  if (outcome == "gaussian") {
    list(alpha = if (base == "sc1") 0.42 else -0.08, delta = 1,
         sigma_rct = 1.5, sigma_rwd = 1.8)
  } else {
    list(alpha = if (base == "sc1") -0.20 else -0.70, delta = log(1.4),
         sigma_rct = 1.4, sigma_rwd = 1.8)
  }
}

# Region indicator for Sc4/Sc5. region = "X5" gives R = {X5 > 2}; "X7" gives
# R = {X7 <= 2}; "none" gives R = empty set (Sc5, shift everywhere).
in_region <- function(X, region) {
  switch(region,
         X5 = X[, 5] > 2,
         X7 = X[, 7] <= 2,
         none = rep(FALSE, nrow(X)),
         stop("unknown region: ", region))
}

# Covariates for a pool of n subjects. Returns X (n x 10), U (n x 2, NA unless
# Sc3) and W (the matrix entering the outcome).
gen_covariates <- function(n, scenario, rho) {
  X <- matrix(rnorm(n * 10), n, 10)
  X <- sweep(X, 2, DGP_MU, "+")
  colnames(X) <- paste0("X", 1:10)
  # 15% of X1 replaced by an extreme value N(+-3, 0.5^2), sign uniform.
  ext <- runif(n) < 0.15
  sgn <- ifelse(runif(n) < 0.5, -1, 1)
  X[ext, 1] <- rnorm(sum(ext), mean = 3 * sgn[ext], sd = 0.5)
  U <- matrix(NA_real_, n, 2, dimnames = list(NULL, c("U1", "U2")))
  W <- X
  if (scenario == 3) {
    U[, 1] <- rnorm(n, DGP_MU[5], 1)
    U[, 2] <- rnorm(n, DGP_MU[6], 1)
    Ustd <- scale(U)
    eps <- matrix(rnorm(n * 2), n, 2)
    X[, 5:6] <- rho * Ustd + sqrt(1 - rho^2) * eps
    W[, 5:6] <- U
  }
  list(X = X, U = U, W = W)
}

# Source indicator D (1 = RCT) for the pool.
gen_source <- function(cov, scenario) {
  n <- nrow(cov$X)
  if (scenario == 1 || scenario >= 4) return(rbinom(n, 1, 0.5))
  V <- if (scenario == 2) cov$X[, 5:6] else cov$U
  ell <- DGP_BETA_D * rowSums(V)
  p <- plogis((ell - median(ell)) / 0.5)
  rbinom(n, 1, p)
}

# Log-normal RMST at horizon tau: E[min(T, tau)] with log T ~ N(mu, sigma^2).
lnorm_rmst <- function(mu, sigma, tau) {
  lt <- log(tau)
  exp(mu + sigma^2 / 2) * pnorm((lt - mu - sigma^2) / sigma) +
    tau * (1 - pnorm((lt - mu) / sigma))
}

# Population median of the standardized log-normal survival mixture over the
# profiles mu (vector) with common sigma: solves mean_i Phi((m - mu_i)/sigma)
# = 0.5 for m on the log scale, returns exp(m).
lnorm_pop_median <- function(mu, sigma) {
  fn <- function(m) mean(pnorm((m - mu) / sigma)) - 0.5
  lo <- min(mu) - 6 * sigma
  hi <- max(mu) + 6 * sigma
  exp(uniroot(fn, c(lo, hi), tol = 1e-8)$root)
}

# Arm-specific truths standardized to a set of profiles (mu0, mu1 vectors).
standardized_truth <- function(mu0, mu1, outcome, sigma, tau = 3) {
  if (outcome == "gaussian") {
    return(list(ate = mean(mu1 - mu0), trt = mean(mu1), ctrl = mean(mu0)))
  }
  med1 <- lnorm_pop_median(mu1, sigma)
  med0 <- lnorm_pop_median(mu0, sigma)
  rm1 <- mean(lnorm_rmst(mu1, sigma, tau))
  rm0 <- mean(lnorm_rmst(mu0, sigma, tau))
  list(ate_med = med1 / med0, trt_med = med1, ctrl_med = med0,
       ate_rmst = rm1 / rm0, trt_rmst = rm1, ctrl_rmst = rm0)
}

#' Generate one replicate.
#'
#' @param scenario integer 1..5. 4 = partial compatibility (RWD shifted by
#'   delta_rwd outside region R), 5 = fully shifted RWD (same as 4 with
#'   region = "none").
#' @param outcome "gaussian" or "survival".
#' @param n_rct,n_rwd RCT pool size and RWD size (exact).
#' @param ratio treatment:control randomization ratio in the RCT (2 for 2:1,
#'   3 for 3:1). Default follows the paper: 2 for Gaussian, 3 for survival.
#' @param rho correlation for Sc3 (ignored otherwise).
#' @param delta_rwd RWD outcome shift for Sc4/Sc5 (Gaussian scale, or log
#'   scale for survival). Ignored for Sc1 to Sc3.
#' @param region "X5" (primary, R = {X5 > 2}), "X7" (secondary, R = {X7 <= 2})
#'   or "none". Forced to "none" for Sc5.
#' @param slope_variant if TRUE, beta5 is multiplied by 0.5 in R^c for the RWD
#'   only (optional Sc4 variant of the spec).
#' @param seed RNG seed (integer). The replicate is a deterministic function
#'   of (all arguments, seed).
#' @param n_mc Monte Carlo pool size for the population truths.
#' @param tau RMST horizon (survival).
#' @return list with data frames rct_trt, rct_ctrl, rwd_ctrl (columns X1..X10,
#'   d, z, and y for Gaussian or time/status for survival), x_rct (all RCT
#'   profiles, treated first then control), truth (sample-standardized, the
#'   primary), truth_pop (population, Monte Carlo), mu0_rct and mu1_rct (true
#'   arm surfaces at the RCT profiles), mu0_rwd_at_rct (the RWD control surface
#'   evaluated at the RCT profiles; equals mu0_rct except in Sc4/Sc5), region_rct
#'   (logical, in R), latent (U and uncensored log times), and constants.
gen_data <- function(scenario, outcome = c("gaussian", "survival"),
                     n_rct = 300, n_rwd = 300, ratio = NULL, rho = 0,
                     delta_rwd = 0, region = "X5", slope_variant = FALSE,
                     seed = 1, n_mc = 1e5, tau = 3) {
  outcome <- match.arg(outcome)
  stopifnot(scenario %in% 1:5)
  if (is.null(ratio)) ratio <- if (outcome == "gaussian") 2 else 3
  if (scenario == 5) region <- "none"
  if (scenario <= 3) { delta_rwd <- 0; slope_variant <- FALSE }
  if (scenario != 3) rho <- 0
  cst <- dgp_constants(scenario, outcome)
  set.seed(seed)

  # Pool, source assignment, exact subsampling.
  n_pool <- 2 * (n_rct + n_rwd)
  repeat {
    cov <- gen_covariates(n_pool, scenario, rho)
    D <- gen_source(cov, scenario)
    if (sum(D == 1) >= n_rct && sum(D == 0) >= n_rwd) break
    n_pool <- 2 * n_pool
  }
  id_rct <- sample(which(D == 1), n_rct)
  id_rwd <- sample(which(D == 0), n_rwd)

  X_rct <- cov$X[id_rct, , drop = FALSE]
  W_rct <- cov$W[id_rct, , drop = FALSE]
  U_rct <- cov$U[id_rct, , drop = FALSE]
  X_rwd <- cov$X[id_rwd, , drop = FALSE]
  W_rwd <- cov$W[id_rwd, , drop = FALSE]
  U_rwd <- cov$U[id_rwd, , drop = FALSE]

  # Treatment assignment in the RCT.
  z <- rbinom(n_rct, 1, ratio / (ratio + 1))

  # Conditional means.
  f_rct <- f_mean(W_rct)
  mu0_rct <- cst$alpha + f_rct
  mu1_rct <- mu0_rct + cst$delta
  f_rwd <- f_mean(W_rwd)
  inR_rwd <- in_region(X_rwd, region)
  inR_rct <- in_region(X_rct, region)
  if (slope_variant) {
    beta_mod <- DGP_BETA; beta_mod[5] <- 0.5 * DGP_BETA[5]
    f_rwd_mod <- f_mean(W_rwd, beta_mod)
    f_rwd <- ifelse(inR_rwd, f_rwd, f_rwd_mod)
    f_rct_mod <- f_mean(W_rct, beta_mod)
  }
  mu_rwd <- cst$alpha + f_rwd + delta_rwd * (!inR_rwd)
  # RWD control surface evaluated at the RCT profiles (for the borrowing map).
  mu0_rwd_at_rct <- cst$alpha +
    (if (slope_variant) ifelse(inR_rct, f_rct, f_rct_mod) else f_rct) +
    delta_rwd * (!inR_rct)

  # Outcomes: the RNG order is fixed so that Sc4 with delta_rwd = 0 equals Sc1.
  y_rct <- mu0_rct + cst$delta * z + rnorm(n_rct, 0, cst$sigma_rct)
  y_rwd <- mu_rwd + rnorm(n_rwd, 0, cst$sigma_rwd)

  make_df <- function(X, y, d, z) {
    df <- as.data.frame(X)
    df$d <- d; df$z <- z
    if (outcome == "gaussian") {
      df$y <- y
    } else {
      n <- length(y)
      t_event <- exp(y)
      c_drop <- rexp(n, rate = 0.02)
      c_admin <- 3 - runif(n)
      c_all <- pmin(c_drop, c_admin)
      df$time <- pmin(t_event, c_all)
      df$status <- as.integer(t_event <= c_all)
    }
    df
  }
  # Order: RCT (treated first, then control), then RWD. Censoring draws are
  # made per data frame in this order.
  ord <- order(-z)
  rct_df <- make_df(X_rct[ord, , drop = FALSE], y_rct[ord], 1L, z[ord])
  rwd_df <- make_df(X_rwd, y_rwd, 0L, 0L)

  # Truths standardized to the replicate's RCT profiles (primary).
  mu0_o <- mu0_rct[ord]; mu1_o <- mu1_rct[ord]
  truth <- standardized_truth(mu0_o, mu1_o, outcome, cst$sigma_rct, tau)

  # Population truths by Monte Carlo over X | D = 1, with an isolated RNG state.
  truth_pop <- NULL
  if (is.finite(n_mc) && n_mc > 0) {
    old <- .Random.seed
    set.seed(seed + 7e6)
    covp <- gen_covariates(n_mc, scenario, rho)
    Dp <- gen_source(covp, scenario)
    Wp <- covp$W[Dp == 1, , drop = FALSE]
    m0 <- cst$alpha + f_mean(Wp)
    truth_pop <- standardized_truth(m0, m0 + cst$delta, outcome, cst$sigma_rct, tau)
    truth_pop$n_rct_pop <- nrow(Wp)
    .Random.seed <<- old
  }

  list(
    rct_trt = rct_df[rct_df$z == 1, ],
    rct_ctrl = rct_df[rct_df$z == 0, ],
    rwd_ctrl = rwd_df,
    x_rct = as.matrix(rct_df[, paste0("X", 1:10)]),
    z_rct = rct_df$z,
    truth = truth,
    truth_pop = truth_pop,
    mu0_rct = mu0_o,
    mu1_rct = mu1_o,
    mu0_rwd_at_rct = mu0_rwd_at_rct[ord],
    region_rct = inR_rct[ord],
    latent = list(U_rct = U_rct[ord, , drop = FALSE], U_rwd = U_rwd,
                  logy_rct = y_rct[ord], logy_rwd = y_rwd),
    constants = c(cst, list(scenario = scenario, outcome = outcome, rho = rho,
                            delta_rwd = delta_rwd, region = region,
                            slope_variant = slope_variant, ratio = ratio,
                            seed = seed, tau = tau))
  )
}

# Named scenario specifications used by the harness. A scenario id string is
# parsed into the gen_data arguments. Examples: "Sc1", "Sc2", "Sc3_rho-0.5",
# "Sc3_rho0", "Sc3_rho0.5", "Sc4_X5_d1", "Sc4_X7_d0.5", "Sc5_d2".
parse_scenario <- function(id) {
  parts <- strsplit(id, "_")[[1]]
  sc <- as.integer(sub("Sc", "", parts[1]))
  out <- list(id = id, scenario = sc, rho = 0, delta_rwd = 0, region = "X5")
  for (p in parts[-1]) {
    if (grepl("^rho", p)) out$rho <- as.numeric(sub("rho", "", p))
    else if (grepl("^d", p)) out$delta_rwd <- as.numeric(sub("d", "", p))
    else if (p %in% c("X5", "X7", "none")) out$region <- p
    else if (p == "slope") out$slope_variant <- TRUE
    else stop("cannot parse scenario token: ", p)
  }
  if (sc == 5) out$region <- "none"
  out
}

gen_data_id <- function(id, outcome, seed, ...) {
  sp <- parse_scenario(id)
  gen_data(scenario = sp$scenario, outcome = outcome, rho = sp$rho,
           delta_rwd = sp$delta_rwd, region = sp$region,
           slope_variant = isTRUE(sp$slope_variant), seed = seed, ...)
}
