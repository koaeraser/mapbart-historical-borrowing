## Trimmed: only ess(mix, sigma=...) for a normMix with the default
## elir method and gaussian/identity family.  Beta/gamma mixtures and
## the morita/moment methods are removed; ess.normMix's binomial /
## poisson family branches are also removed.

#' @export
ess <- function(mix, ...) UseMethod("ess")

#' @export
ess.normMix <- function(mix, sigma, ...) {
  stopifnot(!missing(sigma), is.numeric(sigma), length(sigma) == 1, sigma > 0)
  tauSq <- sigma^2

  ## Gaussian / identity link: Fisher info inverse is 1; we account
  ## for the scale by multiplying the elir integral by sigma^2.
  fisher_inverse <- function(x) 1
  elir <- integrate_density(lir(mix, normMixInfo, fisher_inverse), mix)
  elir * tauSq
}

## Negative second derivative of log mixture density at x, computed
## from the per-component score and curvature in a numerically stable
## log-space form.
mixInfo <- function(mix, x, dens, gradl, hessl) {
  p <- mix[1, ]; a <- mix[2, ]; b <- mix[3, ]
  lp <- log(p)
  ldensComp <- dens(x, a, b, log = TRUE)
  ldensMix  <- matrixStats::logSumExp(lp + ldensComp)
  lwdensComp <- lp + ldensComp - ldensMix
  dgl <- gradl(x, a, b)
  dhl <- hessl(x, a, b) + dgl^2
  if (all(!is.na(dgl)) && (all(dgl < 0) || all(dgl > 0))) {
    gsum <- exp(2 * matrixStats::logSumExp(lwdensComp + log(abs(dgl))))
  } else {
    gsum <- (sum(exp(lwdensComp) * dgl))^2
  }
  if (all(!is.na(dhl)) && (all(dhl < 0) || all(dhl > 0))) {
    hsum <- sign(dhl[1]) *
      exp(matrixStats::logSumExp(lwdensComp + log(abs(dhl))))
  } else {
    hsum <- sum(exp(lwdensComp) * dhl)
  }
  gsum - hsum
}

## Local information ratio: info(x) / fisher_info(x), evaluated
## componentwise; integrate_density takes the expectation under mix.
lir <- function(mix, info, fisher_inverse) {
  Vectorize(function(x) info(mix, x) * fisher_inverse(x))
}

normLogGrad <- function(x, mu, sigma) -(x - mu) / sigma^2
normLogHess <- function(x, mu, sigma) -1 / sigma^2
normMixInfo <- function(mix, x) mixInfo(mix, x, dnorm, normLogGrad, normLogHess)
