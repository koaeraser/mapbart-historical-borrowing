## ============================================================
## helpers.R
##
## Low-level utilities shared across the gMAP -> automixfit -> ess
## pipeline.  Each section is small and self-contained; the entry
## points (gMAP, automixfit, mixfit, ess) and the EM algorithms live
## in their own files and only call into this one.
## ============================================================


## ------------------------------------------------------------
## 1. Functional / array utilities
## ------------------------------------------------------------

## Partial application: Curry(f, x) returns function(...) f(x, ...).
Curry <- function(FUN, ...) {
  .orig <- list(...)
  function(...) do.call(FUN, c(.orig, list(...)))
}

## Return an array shaped like x but filled with `value`, preserving
## class and attributes.  Used to build the trivial identity link.
fill <- function(x, value) {
  ax <- array(value, dim(as.array(x)))
  class(ax) <- class(x)
  attributes(ax) <- attributes(x)
  ax
}


## ------------------------------------------------------------
## 2. Logit / inv_logit (numerically stable variants)
## ------------------------------------------------------------

.bin       <- binomial()
logit      <- .bin$linkfun
inv_logit  <- .bin$linkinv

log_inv_logit <- function(mat) {
  idx <- mat < 0
  mat[idx]  <- mat[idx] - log1p(exp(mat[idx]))
  mat[!idx] <- -log1p(exp(-mat[!idx]))
  mat
}


## ------------------------------------------------------------
## 3. logLik for EM-fitted mixtures (used by automixfit's AIC)
## ------------------------------------------------------------

logLik.EM <- function(object, ...) {
  val <- attr(object, "lli")
  attr(val, "df")   <- attr(object, "df")
  attr(val, "nobs") <- attr(object, "nobs")
  class(val) <- "logLik"
  val
}


## ------------------------------------------------------------
## 4. Density link infrastructure
## ------------------------------------------------------------
##
## A `dlink` is a one-to-one transform attached to a mixture so the
## same density code can serve identity, logit and log links.  We only
## use the identity link here.

dlink <- function(object) attr(object, "link")
"dlink<-" <- function(object, value) {
  attr(object, "link") <- if (is.dlink(value)) value else match.fun(value)()
  object
}

dlink_new <- function(name, link, invlink, lJinv_link) {
  structure(list(name = name, link = link,
                 invlink = invlink, lJinv_link = lJinv_link),
            class = "dlink")
}

is.dlink          <- function(x) inherits(x, "dlink")
is.dlink_identity <- function(x) is.dlink(x) && x$name == "identity"
canonical_dlink   <- function(mix) mix

identity_dlink <- dlink_new("identity",
                            identity, identity,
                            Curry(fill, value = 0))


## ------------------------------------------------------------
## 5. Mixture base: matrix layout, link helpers, [[ subset
## ------------------------------------------------------------

## Build a mixture parameter matrix from per-component vectors and
## tag it with the identity link.  Rows are (weight, p1, p2, ...).
mixdist3 <- function(...) {
  args <- list(...)
  Nc   <- length(args)
  l    <- sapply(args, length)
  if (!all(l >= 3) || !all(l == l[1]))
    stop("All components must have equal number of parameters.")
  res <- do.call(cbind, args)
  if (is.null(names(args)))
    colnames(res) <- paste0("comp", seq(Nc))
  norm <- sum(res[1, ])
  if (norm != 1) res[1, ] <- res[1, ] / norm
  dlink(res) <- identity_dlink
  res
}

mixlink             <- function(mix, x) attr(mix, "link")$link(x)
mixinvlink          <- function(mix, x) attr(mix, "link")$invlink(x)
mixlJinv_link       <- function(mix, l) attr(mix, "link")$lJinv_link(l)
is.mixidentity_link <- function(mix, l) is.dlink_identity(attr(mix, "link"))

## [[ subset a mixture: returns the same class with rescaled weights
## (and preserving the reference scale sigma for normMix).
"[[.mix" <- function(mix, ..., rescale = FALSE) {
  cl <- grep("mix$", class(mix), ignore.case = TRUE, value = TRUE)
  dl <- dlink(mix)
  if (inherits(mix, "normMix")) s <- sigma(mix)
  mix <- mix[, ..., drop = FALSE]
  if (rescale) mix[1, ] <- mix[1, ] / sum(mix[1, ])
  class(mix) <- cl
  dlink(mix) <- dl
  if (inherits(mix, "normMix")) sigma(mix) <- s
  mix
}


## ------------------------------------------------------------
## 6. Attribute setters / accessors on mixtures
## ------------------------------------------------------------

"likelihood<-" <- function(mix, value) {
  attr(mix, "likelihood") <- value
  mix
}

sigma.normMix <- function(object, ...) attr(object, "sigma")

"sigma<-" <- function(object, value) {
  assert_number(value, lower = 0, null.ok = TRUE)
  attr(object, "sigma") <- unname(value)
  object
}


## ------------------------------------------------------------
## 7. Mixture density / CDF / quantile
## ------------------------------------------------------------
##
## Generics dispatch on class; we only register normMix methods.
## The *_impl helpers contain the shared vectorised math.

dmix <- function(mix, x, log = FALSE)                      UseMethod("dmix")
pmix <- function(mix, q, lower.tail = TRUE, log.p = FALSE) UseMethod("pmix")
qmix <- function(mix, p, lower.tail = TRUE, log.p = FALSE) UseMethod("qmix")

dmix_impl <- function(dens, mix, x, log) {
  Nc <- ncol(mix); Nx <- length(x)
  if (is.mixidentity_link(mix)) {
    log_dens <- matrixStats::colLogSumExps(matrix(
      log(mix[1, ]) +
        dens(rep(x, each = Nc),
             rep(mix[2, ], times = Nx),
             rep(mix[3, ], times = Nx),
             log = TRUE),
      nrow = Nc))
  } else {
    ox <- rep(mixinvlink(mix, x), each = Nc)
    log_dens <- matrixStats::colLogSumExps(matrix(
      log(mix[1, ]) +
        rep(mixlJinv_link(mix, x), each = Nc) +
        dens(ox,
             rep(mix[2, ], times = Nx),
             rep(mix[3, ], times = Nx),
             log = TRUE),
      nrow = Nc))
  }
  if (log) log_dens else exp(log_dens)
}

pmix_impl <- function(dist, mix, q, lower.tail = TRUE, log.p = FALSE) {
  Nc <- ncol(mix); Nq <- length(q)
  if (!is.mixidentity_link(mix)) q <- mixinvlink(mix, q)
  comp <- matrix(dist(rep(q, each = Nc),
                      rep(mix[2, ], times = Nq),
                      rep(mix[3, ], times = Nq),
                      lower.tail = lower.tail, log.p = log.p),
                 nrow = Nc)
  if (log.p) matrixStats::colLogSumExps(log(mix[1, ]) + comp)
  else       colSums(mix[1, ] * comp)
}

qmix_impl <- function(quant, mix, p, lower.tail = TRUE, log.p = FALSE) {
  Nc <- ncol(mix)
  if (Nc == 1) {
    return(mixlink(mix, quant(p, mix[2, 1], mix[3, 1],
                              lower.tail = lower.tail, log.p = log.p)))
  }
  ## Bracket the support across components, then bisect via uniroot.
  eps <- 1E-1
  plow  <- if (log.p) min(c(eps, exp(p), 1 - exp(p))) / 2
           else       min(c(eps, p,      1 - p))      / 2
  phigh <- 1 - plow
  qlow  <- mixlink(mix, min(quant(rep.int(plow,  Nc), mix[2, ], mix[3, ])))
  qhigh <- mixlink(mix, max(quant(rep.int(phigh, Nc), mix[2, ], mix[3, ])))
  if (is.infinite(qlow))  qlow  <- -sqrt(.Machine$double.xmax)
  if (is.infinite(qhigh)) qhigh <-  sqrt(.Machine$double.xmax)
  pboundary <- pmix(mix, c(qlow, qhigh), lower.tail = lower.tail, log.p = log.p)
  vapply(seq_along(p), function(i) {
    uniroot(function(x) pmix(mix, x, lower.tail = lower.tail,
                              log.p = log.p) - p[i],
            c(qlow, qhigh),
            f.lower = pboundary[1] - p[i],
            f.upper = pboundary[2] - p[i],
            extendInt = "upX")$root
  }, numeric(1))
}

dmix.normMix <- function(mix, x, log = FALSE) dmix_impl(dnorm, mix, x, log)
pmix.normMix <- function(mix, q, lower.tail = TRUE, log.p = FALSE)
  pmix_impl(pnorm, mix, q, lower.tail, log.p)
qmix.normMix <- function(mix, p, lower.tail = TRUE, log.p = FALSE)
  qmix_impl(qnorm, mix, p, lower.tail, log.p)


## ------------------------------------------------------------
## 8. Multivariate normal mixture constructor
## ------------------------------------------------------------
##
## Used internally by EM_mnmm to package an EM result before EM_nmm
## reclasses it to normMix.  The univariate (Nd=1) case is the only
## one exercised here, but we keep the general shape from upstream.

mixmvnorm <- function(...) {
  mix <- mixdist3(...)
  Nc  <- ncol(mix)
  l   <- nrow(mix)
  p   <- (sqrt(1 + 4 * (l - 1)) - 1) / 2
  assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
  mix <- do.call(mixdist3, lapply(1:Nc, function(co) {
    c(mix[1, co],
      mvnorm(mix[2:(p + 1), co],
             matrix(mix[(p + 2):(1 + p + p^2), co], p, p)))
  }))
  rownames(mix) <- c("w", mvnorm_label(mix[-1, 1], as.character(1:p)))
  class(mix) <- c("mvnormMix", "mix")
  likelihood(mix) <- "mvnormal"
  mix
}

mvnorm <- function(mean, sigma) {
  p <- length(mean)
  s <- sqrt(diag(sigma))
  zero <- s == 0
  rho <- diag(1, nrow = p, ncol = p)
  if (!all(zero))
    rho[!zero, !zero] <- cov2cor(sigma[!zero, !zero, drop = FALSE])
  c(mean, s, rho[lower.tri(rho)])
}

mvnormdim <- function(mvn) as.integer((-3 + sqrt(9 + 8 * length(mvn))) / 2)

mvnorm_label <- function(mvn, dim_labels) {
  p <- mvnormdim(mvn)
  if (p > 1) {
    rho <- outer(dim_labels, dim_labels, paste, sep = ",")[lower.tri(diag(p))]
    c(paste0("m[", dim_labels, "]"),
      paste0("s[", dim_labels, "]"),
      paste0("rho[", rho, "]"))
  } else {
    c(paste0("m[", dim_labels, "]"), paste0("s[", dim_labels, "]"))
  }
}


## ------------------------------------------------------------
## 9. Numerical integration over a mixture density
## ------------------------------------------------------------
##
## Splits the integral by mixture component and substitutes via the
## logit so the integration is on (-Inf, Inf).  Each component's
## quantile function maps back to the natural scale.

integrate_density <- function(integrand, mix,
                              Lplower = -Inf, Lpupper = Inf,
                              eps = getOption("RBesT.integrate_prob_eps", 1E-6)) {
  comp_logit <- function(mc) function(l) {
    u  <- inv_logit(l)
    lp <- log_inv_logit(l)
    ln <- log_inv_logit(-l)
    exp(lp + ln) * integrand(qmix(mc, u))
  }
  Nc <- ncol(mix)
  lower <- inv_logit(Lplower); upper <- inv_logit(Lpupper)
  sum(vapply(1:Nc, function(comp) {
    mc <- mix[[comp, rescale = TRUE]]
    fn <- comp_logit(mc)
    if (all(!is.na(fn(c(Lplower, Lpupper)))))
      return(.integrate(fn, Lplower, Lpupper))
    lo <- if (Lplower == -Inf) qmix(mc, eps)     else qmix(mc, lower)
    hi <- if (Lpupper ==  Inf) qmix(mc, 1 - eps) else qmix(mc, upper)
    .integrate(function(x) integrand(x) * dmix(mc, x), lo, hi)
  }, numeric(1)) * mix[1, ])
}

.integrate <- function(integrand, lower, upper) {
  integrate(integrand, lower = lower, upper = upper,
            rel.tol = .Machine$double.eps^0.25,
            abs.tol = .Machine$double.eps^0.25,
            subdivisions = 1000, stop.on.error = TRUE)$value
}
