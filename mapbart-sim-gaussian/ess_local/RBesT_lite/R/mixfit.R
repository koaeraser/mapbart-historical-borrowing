## Trimmed: only the path used here -- mixfit() generic, the gMAP S3
## method (extracts theta_resp_pred draws from a Stan fit and forwards
## to the default), and the default for type = "norm" (calls EM_nmm).
## Other types (beta, gamma, mvnorm) and the array / gMAPpred methods
## are removed.

#' @export
mixfit <- function(sample, type, thin, ...) UseMethod("mixfit")

#' @export
mixfit.default <- function(sample, type = "norm", thin, ...) {
  stopifnot(identical(type, "norm"))
  if (!missing(thin)) {
    assert_that(thin >= 1)
    sample <- asub(sample, seq(1, NROW(sample), by = thin),
                   dims = 1, drop = FALSE)
  }
  EM_nmm(sample, ...)
}

#' @export
mixfit.gMAP <- function(sample, type, thin, ...) {
  family <- sample$family$family
  if (missing(thin)) thin <- sample$thin
  assert_that(thin >= 1)
  sim <- rstan::extract(sample$fit, pars = "theta_resp_pred",
                        inc_warmup = FALSE, permuted = FALSE)
  sim <- as.vector(sim[seq(1, dim(sim)[1], by = thin), , ])
  mix <- mixfit.default(sim, type = "norm", thin = 1, ...)
  if (!is.null(sample$sigma_ref)) sigma(mix) <- sample$sigma_ref
  set_likelihood(mix, family)
}

set_likelihood <- function(mix, family) {
  if (family == "gaussian") likelihood(mix) <- "normal"
  mix
}
