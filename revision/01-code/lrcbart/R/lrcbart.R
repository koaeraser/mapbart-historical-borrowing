# lrcbart: R interface to the locally robust commensurate BART sampler.
# See method_spec.md (v2) Sections A1, B, C and G, and NOTES.md for choices.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

as_x_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)
  storage.mode(x) <- "double"
  x
}

# Cutpoint grid per variable: numcut equally spaced interior points on
# (min, max) as in BART::wbart with usequants = FALSE; variables with at most
# numcut distinct values use the midpoints between consecutive distinct values.
make_cutpoints <- function(x, numcut = 100) {
  lapply(seq_len(ncol(x)), function(j) {
    u <- sort(unique(x[, j]))
    if (length(u) <= 1) return(numeric(0))
    if (length(u) <= numcut) return((u[-1] + u[-length(u)]) / 2)
    r <- range(u)
    r[1] + (r[2] - r[1]) * seq_len(numcut) / (numcut + 1)
  })
}

# Warm start of the g forest. A g_init is list(nodes, roots): nodes is a
# matrix with columns var (1-based; NA for a leaf), cut (value on the
# covariate scale), left, right (1-based row indices, NA for a leaf), mu, z;
# roots is the 1-based row index of each tree's root (trees beyond the list
# stay stumps). Cut values are snapped to the nearest grid cutpoint.
g_init_to_cpp <- function(g_init, cutpoints) {
  M <- as.matrix(g_init$nodes)
  stopifnot(ncol(M) == 6)
  out <- matrix(-1, nrow(M), 6)
  for (k in seq_len(nrow(M))) {
    if (is.na(M[k, 3]) || M[k, 3] < 1) {
      out[k, ] <- c(-1, -1, -1, -1, M[k, 5], M[k, 6])
    } else {
      v <- M[k, 1]; cp <- cutpoints[[v]]
      if (length(cp) == 0) stop("g_init splits on a variable without cutpoints")
      out[k, ] <- c(v - 1, which.min(abs(cp - M[k, 2])) - 1, M[k, 3] - 1, M[k, 4] - 1, M[k, 5], M[k, 6])
    }
  }
  list(nodes = out, roots = as.integer(g_init$roots) - 1L)
}

#' Greedy regression tree on an offset d(x) (for instance the BART-PP-implied
#' source offset at the RCT control profiles), returned as a g_init for
#' lrc_bart(): one tree of depth at most maxdepth in tree 1, the other H_g - 1
#' trees stumps at zero. Leaves with |mu| > c_slab are marked slab (z = 0).
#' @export
lrc_g_init_greedy <- function(x, d, H_g = 5, maxdepth = 2, n_min = 5, c_slab = 0.5,
                              numcut = 100) {
  x <- as_x_matrix(x); d <- as.numeric(d)
  stopifnot(nrow(x) == length(d))
  cp <- make_cutpoints(x, numcut)
  nodes <- list()
  add <- function(row) { nodes[[length(nodes) + 1]] <<- row; length(nodes) }
  build <- function(idx, depth) {
    mu <- mean(d[idx])
    self <- add(c(NA, NA, NA, NA, mu, as.numeric(abs(mu) <= c_slab)))
    if (depth >= maxdepth || length(idx) < 2 * n_min) return(self)
    best <- NULL; sse0 <- sum((d[idx] - mu)^2)
    for (v in seq_len(ncol(x))) for (c in cp[[v]]) {
      L <- idx[x[idx, v] < c]; R <- idx[x[idx, v] >= c]
      if (length(L) < n_min || length(R) < n_min) next
      sse <- sum((d[L] - mean(d[L]))^2) + sum((d[R] - mean(d[R]))^2)
      if (is.null(best) || sse < best$sse) best <- list(sse = sse, v = v, c = c, L = L, R = R)
    }
    if (is.null(best) || best$sse >= sse0 * 0.999) return(self)
    l <- build(best$L, depth + 1); r <- build(best$R, depth + 1)
    nodes[[self]][1:4] <<- c(best$v, best$c, l, r)
    self
  }
  root <- build(seq_along(d), 0)
  M <- do.call(rbind, nodes)
  colnames(M) <- c("var", "cut", "left", "right", "mu", "z")
  roots <- root
  if (H_g > 1) {
    for (h in 2:H_g) { M <- rbind(M, c(NA, NA, NA, NA, 0, 1)); roots <- c(roots, nrow(M)) }
  }
  list(nodes = M, roots = roots)
}

# Residual SD guess for the Inv-Gamma scale (paper: OLS residual SD, or a
# parametric log-normal AFT fit when there is censoring).
sigma_hat <- function(y, x, status = NULL) {
  n <- length(y)
  if (n < 2) return(1)
  if (!is.null(status) && any(status == 0)) {
    if (requireNamespace("survival", quietly = TRUE) && n > ncol(x) + 5) {
      fit <- tryCatch(
        survival::survreg(survival::Surv(exp(y), status) ~ x, dist = "lognormal"),
        error = function(e) NULL)
      if (!is.null(fit) && is.finite(fit$scale)) return(fit$scale)
    }
    return(stats::sd(y[status == 1]))
  }
  if (n > ncol(x) + 5) {
    fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)$sigma
      if (is.finite(s) && s > 0) return(s)
    }
  }
  stats::sd(y)
}

# ---------------------------------------------------------------------------
# core fitting routine (shared by lrc_bart, lrc_bart_aft, bart_fit)
# ---------------------------------------------------------------------------

lrc_core <- function(y, x, source, status = NULL, aft = FALSE,
                     H_f = 50, H_g = 5,
                     alpha_f = 0.95, beta_f = 2, alpha_g = 0.5, beta_g = 3,
                     k = 2, nu_sigma = 3, sigquant = 0.9,
                     nu0 = 3, s0_sq = 0.01, update_tau0 = TRUE,
                     w = NULL, a_w = 1, b_w = 1,
                     tau1_sq = NULL, tau1_mult = 1, update_tau1 = FALSE, nu1 = 3, s1_sq = NULL,
                     n_min = 5, n_min_g = n_min, numcut = 100,
                     p_grow = 0.3, p_prune = 0.3,
                     n_burn = 1000, n_draw = 1000, thin = 1, g_sweeps = 1, g_init = NULL,
                     seed = NULL, keep_trees = TRUE, keep_train = TRUE,
                     keep_leaf_stats = FALSE,
                     verbose = FALSE, calib = NULL, single_arm = NULL) {
  x <- as_x_matrix(x)
  n <- nrow(x)
  y <- as.numeric(y)
  stopifnot(length(y) == n, length(source) == n)
  if (is.character(source) || is.factor(source)) {
    source <- ifelse(tolower(as.character(source)) %in% c("rct", "1", "rct_control"), 1L, 2L)
  }
  source <- as.integer(source)
  stopifnot(all(source %in% c(1L, 2L)))
  if (is.null(status)) status <- rep(1L, n)
  status <- as.integer(status)
  stopifnot(all(status %in% c(0L, 1L)))
  if (!aft && any(status == 0)) stop("censoring requires lrc_bart_aft()")
  if (any(!is.finite(y))) stop("y must be finite")
  n1 <- sum(source == 1L)
  is_single <- n1 == 0
  if (!is.null(single_arm) && single_arm != is_single)
    stop("single_arm = ", single_arm, " but the data have ", n1, " RCT controls")
  if (!is.null(seed)) set.seed(seed)

  # centering and range rules on the observed values (log times for AFT)
  y_center <- mean(y)
  yc <- y - y_center
  rng <- diff(range(yc))
  lambda_f_sq <- (rng / (2 * k * sqrt(max(H_f, 1))))^2
  # slab scale: the range rule for an H_g-tree ensemble times tau1_mult
  if (is.null(tau1_sq)) tau1_sq <- tau1_mult * rng^2 / (4 * k^2 * max(H_g, 1))
  if (is.null(s1_sq)) s1_sq <- tau1_sq
  if (!is.null(w)) {
    stopifnot(length(w) == 1, w >= 0, w <= 1)
    w_fixed <- w
  } else {
    w_fixed <- NA_real_
  }
  if (!is.null(calib) && is.null(s0_sq)) s0_sq <- calib$s0_sq
  if (is.null(s0_sq)) stop("s0_sq must be given (or a calib object)")

  # residual-variance priors per source
  sh2 <- if (sum(source == 2L) > 0) sigma_hat(yc[source == 2L], x[source == 2L, , drop = FALSE], status[source == 2L]) else NA
  sh1 <- if (n1 > 0) sigma_hat(yc[source == 1L], x[source == 1L, , drop = FALSE], status[source == 1L]) else sh2
  if (is.na(sh2)) sh2 <- sh1
  qchi <- stats::qchisq(1 - sigquant, nu_sigma) / nu_sigma
  lambda_sigma1 <- sh1^2 * qchi
  lambda_sigma2 <- sh2^2 * qchi

  cutpoints <- make_cutpoints(x, numcut)
  logC <- ifelse(status == 0L, yc, -Inf)

  hyper <- list(H_f = as.integer(H_f), H_g = as.integer(H_g),
                lambda_f_sq = lambda_f_sq, nu_sigma = nu_sigma,
                lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2,
                nu0 = nu0, s0_sq = s0_sq, nu1 = nu1, s1_sq = s1_sq,
                a_w = a_w, b_w = b_w, tau1_sq = tau1_sq,
                update_tau0 = update_tau0, update_tau1 = update_tau1,
                aft = aft, w_fixed = w_fixed,
                sigma1_sq_init = sh1^2, sigma2_sq_init = sh2^2,
                p_grow = p_grow, p_prune = p_prune, n_min = as.integer(n_min),
                n_min_g = as.integer(n_min_g),
                alpha_f = alpha_f, beta_f = beta_f, alpha_g = alpha_g, beta_g = beta_g,
                tau1_mult = tau1_mult)
  control <- list(n_burn = as.integer(n_burn), n_draw = as.integer(n_draw),
                  thin = as.integer(thin), verbose = verbose,
                  keep_trees = keep_trees, keep_train = keep_train,
                  keep_leaf_stats = keep_leaf_stats, g_sweeps = as.integer(g_sweeps),
                  g_init = if (is.null(g_init)) NULL else g_init_to_cpp(g_init, cutpoints))

  t0 <- proc.time()[["elapsed"]]
  out <- lrc_sampler_cpp(x, yc, source, status, logC, cutpoints, hyper, control)
  elapsed <- proc.time()[["elapsed"]] - t0

  fit <- list(
    draws = list(f = out$f_draws, g = out$g_draws),
    sigma1_sq = out$sigma1_sq, sigma2_sq = out$sigma2_sq,
    tau0_sq = out$tau0_sq, tau1_sq = out$tau1_sq, w = out$w,
    K0 = out$K0, L_g = out$L_g,
    f_train = if (keep_train) out$f_train + y_center else NULL,
    g_train = if (keep_train && H_g > 0) out$g_train else NULL,
    accept = out$accept,
    g_leaf_stats = if (keep_leaf_stats && H_g > 0) lapply(out$g_leaf_stats, function(m) {
      colnames(m) <- c("tree", "n1", "ybar", "v1", "z", "theta", "pspike"); m }) else NULL,
    y_center = y_center, H_f = H_f, H_g = H_g, aft = aft,
    single_arm = out$single_arm, n1 = n1, n2 = n - n1,
    hyper = hyper, calib = calib, x_train = x, source = source,
    y = y, status = status, n_draw = n_draw, time = elapsed)
  dimnames(fit$accept) <- list(c("f", "g"), c("grow", "prune", "change"))
  class(fit) <- "lrcbart"
  fit
}

#' Locally robust commensurate BART, Gaussian outcome.
#'
#' @param y numeric outcome (RCT controls and RWD controls pooled).
#' @param x covariate matrix or data frame (numeric columns).
#' @param source 1 (or "rct") for RCT controls, 2 (or "rwd") for RWD controls.
#'   No RCT controls at all gives the single-arm mode of spec Section A1.
#' @param w NULL to learn w ~ Beta(a_w, b_w); a number in [0, 1] fixes it
#'   (w = 1 gives C-BART, w = 0 the free-offset model).
#' @param s0_sq scale of the scaled-Inv-chi2(nu0, s0_sq) prior on tau0^2,
#'   normally from lrc_ess_calibrate(); or pass calib = that object.
#' @param ... see lrc_core for the remaining hyperparameters.
#' @export
lrc_bart <- function(y, x, source, ...) {
  lrc_core(y = y, x = x, source = source, status = NULL, aft = FALSE, ...)
}

#' LRC-BART for censored log-normal (AFT) outcomes.
#'
#' @param time positive event or censoring times; @param status 1 event, 0 censored.
#' @export
lrc_bart_aft <- function(time, status, x, source, ...) {
  stopifnot(all(time > 0))
  lrc_core(y = log(time), x = x, source = source, status = status, aft = TRUE, ...)
}

#' Plain BART (the treatment arm, or Stage 1 of the ESS calibration): the same
#' f code path with one source and no g ensemble.
#' @export
bart_fit <- function(y, x, status = NULL, time = NULL, ...) {
  aft <- !is.null(time)
  if (aft) { stopifnot(all(time > 0)); y <- log(time) }
  n <- if (is.matrix(x) || is.data.frame(x)) nrow(x) else length(x)
  fit <- lrc_core(y = y, x = x, source = rep(1L, n), status = status, aft = aft,
                  H_g = 0, s0_sq = 1, ...)
  fit$plain_bart <- TRUE
  fit
}

# ---------------------------------------------------------------------------
# prediction
# ---------------------------------------------------------------------------

#' Posterior draws of the control surface at new covariates.
#' @param arm "rct_control" gives f + g, "rwd_control" gives f.
#' @return matrix n_draw x nrow(x_new).
#' @export
predict.lrcbart <- function(object, x_new, arm = c("rct_control", "rwd_control"), ...) {
  arm <- match.arg(arm)
  x_new <- as_x_matrix(x_new)
  if (is.null(object$draws$f) || length(object$draws$f) == 0)
    stop("fit was run with keep_trees = FALSE")
  f <- predict_forest_cpp(object$draws$f, x_new)$value + object$y_center
  if (arm == "rwd_control" || object$H_g == 0) return(f)
  f + predict_forest_cpp(object$draws$g, x_new)$value
}

#' Posterior draws of the discrepancy g at new covariates.
#' @export
lrc_discrepancy <- function(fit, x_new) {
  if (fit$H_g == 0) stop("no g ensemble in this fit")
  predict_forest_cpp(fit$draws$g, as_x_matrix(x_new))$value
}

#' Borrowing (compatibility) map at new covariates.
#'
#' type = "abs_g" (default): B(x) = P(|g(x)| < c | data), the posterior
#' probability that the outcome-scale discrepancy at x is within c; c is on
#' the outcome scale for Gaussian fits and on the log-time scale for AFT fits
#' (default 0.5, one third of the RCT residual SD in the Gaussian study).
#' type = "all_spike": the ensemble-level map P(all H_g leaves covering x are
#' spike | data), kept as a diagnostic (its scale depends on H_g).
#' @return list(draws = n_draw x n_new 0/1 matrix, mean = posterior mean,
#'   type, c; for "all_spike" also nspike = number of spike trees per draw and point).
#' @export
lrc_borrowing_map <- function(fit, x_new, type = c("abs_g", "all_spike"), c = 0.5) {
  type <- match.arg(type)
  if (fit$H_g == 0) stop("no g ensemble in this fit")
  x_new <- as_x_matrix(x_new)
  if (type == "abs_g") {
    stopifnot(length(c) == 1, is.finite(c), c > 0)
    g <- predict_forest_cpp(fit$draws$g, x_new)$value
    d <- (abs(g) < c) * 1
    return(list(draws = d, mean = colMeans(d), type = type, c = c))
  }
  ns <- predict_forest_cpp(fit$draws$g, x_new)$nspike
  d <- (ns == fit$H_g) * 1
  list(draws = d, mean = colMeans(d), nspike = ns, type = type, c = NA_real_)
}

# ---------------------------------------------------------------------------
# ESS calibration (spec Section C)
# ---------------------------------------------------------------------------

#' Monte Carlo c_g = E[N^-2 sum_{i,j} 1{x_i, x_j share a leaf}] under the
#' g-tree prior, on the cutpoint grid built from x_grid (defaults to x_rct).
#' @export
lrc_cg <- function(x_rct, x_grid = NULL, alpha_g = 0.5, beta_g = 3, n_sim = 2000,
                   numcut = 100, maxdepth = 12, seed = NULL) {
  x_rct <- as_x_matrix(x_rct)
  if (is.null(x_grid)) x_grid <- x_rct
  if (!is.null(seed)) set.seed(seed)
  cg_montecarlo_cpp(x_rct, make_cutpoints(as_x_matrix(x_grid), numcut),
                    alpha_g, beta_g, as.integer(n_sim), as.integer(maxdepth))
}

#' ESS_0(s0_sq) for one block: sigma1^2 E_{tau0^2}[1 / (V + c_g H_g tau0^2)].
ess0_block <- function(V, sigma1_sq, c_g, H_g, tau0_draws) {
  sigma1_sq * mean(1 / (V + c_g * H_g * tau0_draws))
}

#' Prior ESS calibration by the Pr-rule.
#'
#' Stage 1: standard BART on the RWD alone, evaluated at the RCT profiles.
#' V_mu^f = Var(mean_i f(x_i)) is computed from the whole chain; the Pr-rule
#' blocks (n_blocks consecutive blocks of draws) carry only the block-to-block
#' spread, rescaled so that their mean equals the whole-chain variance
#' (within-block variances are biased low by the chain's autocorrelation).
#' Stage 2: ESS_0(s0_sq) on a grid, G(s0_sq) = fraction of blocks with
#' ESS_0 <= N_target, and s0_sq* = smallest grid value with G >= q.
#'
#' sigma_1^2 (the unit Fisher information of an RCT control) is, in order of
#' preference: the argument sigma1_sq; the residual variance of the RCT
#' controls (y_rct_ctrl, x_rct_ctrl, status_rct_ctrl), by default the
#' posterior mean sigma^2 of a plain BART fit on the RCT controls alone
#' (sigma1_method = "bart"; "residual" uses the OLS or log-normal survreg
#' rule of the sampler's residual prior, which keeps the curvature of f in
#' the residual); and, when neither is given (single-arm setting, no RCT
#' outcomes yet), the Stage-1 posterior mean of the RWD residual variance.
#'
#' Feasibility: ESS_0 cannot exceed the ceiling sigma_1^2 / V_mu^f (T3). If
#' N_target is at or above it, N_target is capped at the ceiling, s0_sq is the
#' smallest grid value (the point-mass limit) and ess_capped = TRUE is
#' returned; one warning is issued per calibration when warn = TRUE.
#'
#' @param y_rwd,x_rwd RWD outcome (log time for AFT, with status_rwd) and covariates.
#' @param x_rct the N RCT covariate profiles the estimand is standardized to.
#' @param sigma1_sq RCT residual variance (overrides y_rct_ctrl).
#' @param y_rct_ctrl,x_rct_ctrl,status_rct_ctrl the RCT controls (outcome on
#'   the same scale as y_rwd), used to estimate sigma1_sq when it is not given.
#' @param fit_stage1 an existing bart_fit() on the RWD (skips Stage 1).
#' @export
lrc_ess_calibrate <- function(y_rwd, x_rwd, x_rct, N_target, grid = NULL, q = 0.95,
                              status_rwd = NULL, sigma1_sq = NULL,
                              y_rct_ctrl = NULL, x_rct_ctrl = NULL, status_rct_ctrl = NULL,
                              nu0 = 3, H_g = 5, alpha_g = 0.5, beta_g = 3,
                              H_f = 50, n_blocks = 20, n_burn = 1000, n_draw = 1000,
                              n_tau = 4000, n_cg = 2000, seed = 1,
                              fit_stage1 = NULL, warn = TRUE,
                              sigma1_method = c("bart", "residual"), ...) {
  x_rwd <- as_x_matrix(x_rwd); x_rct <- as_x_matrix(x_rct)
  if (is.null(grid)) grid <- 10^seq(-6, 1, length.out = 141)
  grid <- sort(grid)
  set.seed(seed)
  if (is.null(fit_stage1)) {
    fit_stage1 <- if (is.null(status_rwd) || all(status_rwd == 1))
      bart_fit(y_rwd, x_rwd, H_f = H_f, n_burn = n_burn, n_draw = n_draw, ...)
    else
      bart_fit(y_rwd, x_rwd, status = status_rwd, time = exp(y_rwd), H_f = H_f,
               n_burn = n_burn, n_draw = n_draw, ...)
  }
  # sigma_1^2: RCT residual variance when RCT controls exist, else the RWD's
  sigma1_method <- match.arg(sigma1_method)
  sigma1_source <- "given"
  if (is.null(sigma1_sq)) {
    if (!is.null(y_rct_ctrl) && length(y_rct_ctrl) > 0) {
      x1 <- as_x_matrix(x_rct_ctrl)
      stopifnot(nrow(x1) == length(y_rct_ctrl))
      st1 <- if (is.null(status_rct_ctrl)) NULL else as.integer(status_rct_ctrl)
      if (sigma1_method == "bart") {
        # plain BART on the RCT controls alone: the residual variance of a
        # nonparametric fit (an OLS residual keeps the curvature of f)
        b1 <- if (is.null(st1) || all(st1 == 1))
          bart_fit(y_rct_ctrl, x1, H_f = H_f, n_burn = n_burn, n_draw = n_draw, seed = seed + 11L, keep_trees = FALSE, keep_train = FALSE)
        else
          bart_fit(y_rct_ctrl, x1, status = st1, time = exp(y_rct_ctrl), H_f = H_f, n_burn = n_burn, n_draw = n_draw, seed = seed + 11L, keep_trees = FALSE, keep_train = FALSE)
        sigma1_sq <- mean(b1$sigma1_sq)
        sigma1_source <- "rct_bart"
      } else {
        sigma1_sq <- sigma_hat(y_rct_ctrl - mean(y_rct_ctrl), x1, st1)^2
        sigma1_source <- "rct_residual"
      }
    } else {
      sigma1_sq <- mean(fit_stage1$sigma1_sq)  # the single source of bart_fit is the RWD
      sigma1_source <- "rwd_posterior"
    }
  }
  fdraws <- predict(fit_stage1, x_rct)
  mu_draws <- rowMeans(fdraws)
  nd <- length(mu_draws)
  V_whole <- stats::var(mu_draws)
  blk <- ceiling(seq_len(nd) / (nd / n_blocks))
  V_b_raw <- tapply(mu_draws, blk, stats::var)
  V_b <- V_b_raw * V_whole / mean(V_b_raw)  # block spread around the whole-chain level
  sig1_b <- rep(sigma1_sq, length(V_b))
  c_g <- lrc_cg(x_rct, x_grid = rbind(x_rwd, x_rct), alpha_g = alpha_g, beta_g = beta_g,
                n_sim = n_cg)
  ess_ceiling <- sigma1_sq / V_whole
  ess_capped <- FALSE
  # N_target: an absolute ESS, or c("frac", p) / list("frac", p) for the
  # fraction p of the ceiling sigma1^2 / V_mu^f (a target that is always
  # feasible when p < 1)
  N_target_spec <- N_target
  if (is.list(N_target) || is.character(N_target)) {
    stopifnot(length(N_target) == 2, identical(as.character(N_target[[1]]), "frac"))
    frac <- as.numeric(N_target[[2]])
    stopifnot(is.finite(frac), frac > 0)
    N_target <- frac * ess_ceiling
  }
  N_target_req <- N_target
  if (N_target >= ess_ceiling) {
    ess_capped <- TRUE
    N_target <- ess_ceiling
    if (warn)
      warning(sprintf("N_target %.1f is at or above the ESS ceiling sigma1^2 / V_mu^f = %.1f; capped at the ceiling (s0_sq at the point-mass limit)",
                      N_target_req, ess_ceiling), call. = FALSE)
  }
  chi <- stats::rchisq(n_tau, nu0)  # common random numbers across the grid
  ess_mat <- sapply(grid, function(s0) {
    tau0 <- nu0 * s0 / chi
    mapply(function(V, s1) ess0_block(V, s1, c_g, H_g, tau0), V_b, sig1_b)
  })
  ess_mat <- matrix(ess_mat, nrow = length(V_b))
  G <- colMeans(ess_mat <= N_target)
  ceiling_b <- sig1_b / V_b
  at_boundary <- NULL
  if (ess_capped) {
    s0_sq <- grid[1]; at_boundary <- "lower"
  } else {
    ok <- which(G >= q)
    if (length(ok) == 0) {
      s0_sq <- grid[length(grid)]; at_boundary <- "upper"
      if (warn) warning("no grid value reaches q; returning the largest grid value", call. = FALSE)
    } else {
      s0_sq <- grid[ok[1]]
      if (ok[1] == 1) at_boundary <- "lower"
    }
  }
  list(s0_sq = s0_sq, N_target = N_target, N_target_requested = N_target_req,
       N_target_spec = N_target_spec,
       ess_capped = ess_capped, q = q, grid = grid, G = G,
       ess0_grid = colMeans(ess_mat), ess0_blocks = ess_mat,
       ess0_at_s0 = mean(ess_mat[, match(s0_sq, grid)]),
       ceiling = ess_ceiling, ceiling_blocks = ceiling_b,
       V_mu_f = V_whole, V_blocks = V_b, V_blocks_raw = V_b_raw,
       sigma1_sq = sigma1_sq, sigma1_source = sigma1_source, sigma1_sq_blocks = sig1_b,
       c_g = c_g, H_g = H_g, nu0 = nu0, at_boundary = at_boundary,
       fit_stage1 = fit_stage1)
}

#' Prior-averaged ESS(w, s0_sq) of spec Section C (binomial mixture over the
#' number of spike trees), for reporting alongside ESS_0.
#' @export
lrc_ess_prior <- function(calib, s0_sq = calib$s0_sq, w = 0.5, tau1_sq, n_tau = 4000, seed = 1) {
  set.seed(seed)
  chi <- stats::rchisq(n_tau, calib$nu0)
  tau0 <- calib$nu0 * s0_sq / chi
  H <- calib$H_g
  vals <- sapply(0:H, function(kk) {
    mean(mapply(function(V, s1) s1 * mean(1 / (V + calib$c_g * (kk * tau0 + (H - kk) * tau1_sq))),
                calib$V_blocks, calib$sigma1_sq_blocks))
  })
  sum(stats::dbinom(0:H, H, w) * vals)
}

#' Realized ESS: E_post[sigma1^2 / {V_mu^f + c_g sum_h tau^2_{z_h(x)}}] with
#' (z, tau0^2) from the posterior, averaged over the RCT profiles x_rct.
#' V_mu^f and c_g come from the calibration object stored in the fit (or given);
#' otherwise V_mu^f is approximated by the variance of the fit's own f draws at
#' x_rct (a joint-fit quantity, flagged in the output).
#' @export
lrc_ess_realized <- function(fit, x_rct = NULL, V_mu_f = NULL, c_g = NULL) {
  if (fit$H_g == 0) stop("no g ensemble in this fit")
  if (is.null(x_rct)) x_rct <- fit$x_train[fit$source == 1L, , drop = FALSE]
  x_rct <- as_x_matrix(x_rct)
  approx <- FALSE
  if (is.null(V_mu_f)) {
    if (!is.null(fit$calib)) V_mu_f <- fit$calib$V_mu_f
    else { V_mu_f <- stats::var(rowMeans(predict(fit, x_rct, arm = "rwd_control"))); approx <- TRUE }
  }
  if (is.null(c_g)) {
    c_g <- if (!is.null(fit$calib)) fit$calib$c_g else lrc_cg(x_rct, alpha_g = fit$hyper$alpha_g, beta_g = fit$hyper$beta_g)
  }
  # exact estimand-level discrepancy variance given the posterior partition
  # and indicators: sum_h sum_l (n_l / N)^2 tau^2_{z_hl}; c_g is only used
  # for the spec's tree-level approximation, reported alongside
  Vg <- forest_gvar_cpp(fit$draws$g, x_rct, fit$tau0_sq, fit$tau1_sq)
  draws <- fit$sigma1_sq / (V_mu_f + Vg)
  ns <- predict_forest_cpp(fit$draws$g, x_rct)$nspike
  tbar <- rowMeans(ns * fit$tau0_sq + (fit$H_g - ns) * fit$tau1_sq)
  draws_cg <- fit$sigma1_sq / (V_mu_f + c_g * tbar)
  list(draws = draws, mean = mean(draws), V_g = Vg, V_mu_f = V_mu_f, c_g = c_g,
       mean_cg_approx = mean(draws_cg), V_mu_f_approx = approx)
}

#' Pointwise ESS map ESS(x) = sigma1^2 / {V_f(x) + H_g tau0^2}, with V_f(x) the
#' RWD posterior variance at x (spec Section C, last paragraph).
#' @export
lrc_ess_map <- function(fit_stage1, x_new, tau0_sq, H_g = 5, sigma1_sq = NULL) {
  f <- predict(fit_stage1, x_new)
  Vf <- apply(f, 2, stats::var)
  s1 <- if (is.null(sigma1_sq)) mean(fit_stage1$sigma1_sq) else sigma1_sq
  s1 / (Vf + H_g * mean(tau0_sq))
}

# ---------------------------------------------------------------------------
# standardized ATE (Gaussian difference; survival median and RMST ratios)
# ---------------------------------------------------------------------------

# closed-form log-normal RMST (paper eq. rmst-closed-form), vectorized over a
# matrix mu (draws x profiles) and a vector sigma (per draw)
lnorm_rmst_mat <- function(mu, sigma, tau) {
  lt <- log(tau)
  s <- matrix(sigma, nrow(mu), ncol(mu))
  exp(mu + s^2 / 2) * stats::pnorm((lt - mu - s^2) / s) + tau * stats::pnorm((mu - lt) / s)
}

# population median per draw: solve mean_i Phi((m - mu_i)/sigma) = 0.5 by
# bisection, vectorized over draws
lnorm_pop_median_mat <- function(mu, sigma, iter = 60) {
  lo <- apply(mu, 1, min) - 8 * sigma
  hi <- apply(mu, 1, max) + 8 * sigma
  for (it in seq_len(iter)) {
    mid <- (lo + hi) / 2
    val <- rowMeans(stats::pnorm((mid - mu) / sigma))
    up <- val < 0.5
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  exp((lo + hi) / 2)
}

#' Standardized ATE over the RCT profiles x_rct.
#' @param fit_trt bart_fit() on the treated arm; fit_ctrl lrc_bart() or
#'   lrc_bart_aft() (or bart_fit()) on the controls; the RCT control surface
#'   f + g is used.
#' @return Gaussian: list(ate, trt, ctrl) draws. Survival: list(ate_med,
#'   trt_med, ctrl_med, ate_rmst, trt_rmst, ctrl_rmst) draws.
#' @export
lrc_ate <- function(fit_trt, fit_ctrl, x_rct, outcome = c("gaussian", "survival"), tau = 3) {
  outcome <- match.arg(outcome)
  x_rct <- as_x_matrix(x_rct)
  mu1 <- predict(fit_trt, x_rct)
  mu0 <- predict(fit_ctrl, x_rct, arm = "rct_control")
  nd <- min(nrow(mu1), nrow(mu0))
  mu1 <- mu1[seq_len(nd), , drop = FALSE]; mu0 <- mu0[seq_len(nd), , drop = FALSE]
  if (outcome == "gaussian") {
    trt <- rowMeans(mu1); ctrl <- rowMeans(mu0)
    return(list(ate = trt - ctrl, trt = trt, ctrl = ctrl))
  }
  s1 <- sqrt(fit_trt$sigma1_sq[seq_len(nd)])
  s0 <- sqrt(fit_ctrl$sigma1_sq[seq_len(nd)])
  med1 <- lnorm_pop_median_mat(mu1, s1); med0 <- lnorm_pop_median_mat(mu0, s0)
  rm1 <- rowMeans(lnorm_rmst_mat(mu1, s1, tau)); rm0 <- rowMeans(lnorm_rmst_mat(mu0, s0, tau))
  list(ate_med = med1 / med0, trt_med = med1, ctrl_med = med0,
       ate_rmst = rm1 / rm0, trt_rmst = rm1, ctrl_rmst = rm0)
}

#' @export
print.lrcbart <- function(x, ...) {
  cat("lrcbart fit:", if (isTRUE(x$plain_bart)) "plain BART" else if (x$single_arm) "single-arm LRC-BART" else "two-arm LRC-BART", "\n")
  cat("  n_rct_ctrl =", x$n1, " n_rwd =", x$n2, " H_f =", x$H_f, " H_g =", x$H_g,
      " aft =", x$aft, " draws =", x$n_draw, "\n")
  cat("  sigma1 =", round(sqrt(mean(x$sigma1_sq)), 3), " sigma2 =", round(sqrt(mean(x$sigma2_sq)), 3), "\n")
  if (x$H_g > 0) cat("  tau0_sq =", signif(mean(x$tau0_sq), 3), " tau1_sq =", signif(mean(x$tau1_sq), 3),
                     " w =", round(mean(x$w), 3), " spike fraction =", round(mean(x$K0 / pmax(x$L_g, 1)), 3), "\n")
  cat("  time =", round(x$time, 1), "s\n")
  invisible(x)
}
