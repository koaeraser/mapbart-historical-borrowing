## ============================================================
## EM.R
##
## EM fit of a univariate normal mixture.  The flow is:
##
##   knn       -- k-means clustering, seeds the t-mixture EM
##   EM_msmm   -- multivariate Student-t mixture EM, robust init
##   EM_mnmm   -- multivariate normal mixture EM (uses EM_msmm init)
##   EM_nmm    -- univariate wrapper called by mixfit(type = "norm")
##
## EM_nmm is the only function called from outside this file.  The
## others are kept here because they form a single self-contained
## algorithm and are not used elsewhere.
## ============================================================


## ------------------------------------------------------------
## knn -- k-means clustering used by EM_msmm to seed the EM init
## ------------------------------------------------------------

knn <- function(X, K = 2, init, Ninit = 50, verbose = FALSE,
                tol, Niter.max = 500) {
  if (!is.matrix(X)) X <- matrix(X, ncol = 1)
  if (missing(tol)) tol <- min(colVars(X)) / 5
  N  <- dim(X)[1]
  Nd <- dim(X)[2]

  if (missing(init)) {
    pEst  <- runif(K) / K; pEst <- pEst / sum(pEst)
    muEst <- matrix(0, K, Nd)
    for (i in seq(K))
      muEst[i, ] <- colMeans(X[sample.int(N, min(Ninit, N), replace = FALSE),
                                , drop = FALSE])
  } else {
    pEst  <- init$p
    muEst <- init$mu
  }

  DM <- matrix(0, N, K)
  J  <- Inf
  iter <- 1
  while (iter < Niter.max) {
    Jprev <- J
    for (i in seq(K))
      DM[, i] <- rowSums(scale(X, muEst[i, ], FALSE)^2)
    resp  <- 1 * (1 == matrixStats::rowRanks(DM, ties.method = "first"))
    respM <- matrixStats::colSums2(resp)
    if (any(respM == 0)) {
      warning("Some components assigned empty set; try reducing K.")
      respM[respM == 0] <- 1
    }
    muEst <- sweep(crossprod(resp, X), 1, respM, "/", FALSE)
    J     <- sum((X - resp %*% muEst)^2)
    if (Jprev - J < tol) break
    iter <- iter + 1
  }

  list(center = muEst,
       p = colMeans(resp),
       cluster = ((which(t(resp == 1)) - 1) %% K) + 1)
}


## ------------------------------------------------------------
## EM_msmm -- Student-t mixture EM (robust init for EM_mnmm).
## Reference: Peel & McLachlan (2000).  Caller only reads
## $p, $center, $cov from the returned list.
## ------------------------------------------------------------

EM_msmm <- function(X, Nc, init, Ninit = 50, verbose = FALSE,
                    Niter.max = 500, tol = 1e-1) {

  if (!is.matrix(X)) X <- matrix(X, ncol = 1)
  N  <- dim(X)[1]
  Nd <- dim(X)[2]

  ## --- initialisation: stratified k-means seeds, then knn ---------
  if (missing(init)) {
    ind <- seq(1, N - Nc, length = Ninit)
    knnInit <- list(mu = matrix(0, nrow = Nc, ncol = Nd),
                    p  = (1 / seq(1, Nc)) / sum(seq(1, Nc)))
    for (k in seq(Nc))
      knnInit$mu[k, ] <- colMeans(X[ind + k - 1, , drop = FALSE])
    suppressWarnings(KNN <- knn(X, K = Nc, init = knnInit,
                                verbose = verbose, Niter.max = 50))
    pEst   <- KNN$p
    cmin   <- which.min(pEst)
    muEst  <- KNN$center
    nuEst  <- rep(12, times = Nc)
    covEst <- array(0, dim = c(Nc, Nd, Nd))
    Xtau   <- sqrt(colVars(X))
    for (i in seq(Nc)) {
      if (i == cmin) next
      ind <- KNN$cluster == i
      if (sum(ind) > 10) {
        covKNN <- as.matrix(cov(X[ind, , drop = FALSE]))
        R   <- cov2cor(covKNN)
        tau <- sqrt(diag(covKNN))
        tau[tau <= 0] <- Xtau[tau <= 0]
      } else {
        R   <- diag(Nd)
        tau <- Xtau
      }
      tau <- pmax(tau, Xtau / 100)
      covEst[i, , ] <- diag(tau, Nd, Nd) %*% R %*% diag(tau, Nd, Nd)
    }
    covEst[cmin, , ] <- diag(Xtau, Nd, Nd) %*% diag(Nd) %*% diag(Xtau, Nd, Nd)
    muEst[cmin] <- sum(pEst * muEst)
  } else {
    pEst   <- init$p
    muEst  <- init$center
    nuEst  <- init$nu
    covEst <- init$cov
  }

  iter     <- 0
  logN     <- log(N)
  traceLli <- c()
  Dlli     <- Inf
  lli      <- array(-20, dim = c(N, Nc))
  lU       <- array(0,   dim = c(N, Nc))

  ## degrees-of-freedom MLE objective per component
  nu_ml <- function(c1) function(nu) {
    (log(nu / 2) - digamma(nu / 2) + c1 +
       digamma((nu + Nd) / 2) - log((nu + Nd) / 2))^2
  }

  while (iter < Niter.max) {
    for (i in seq(Nc)) {
      lli[, i] <- log(pEst[i]) +
        mvtnorm::dmvt(X, muEst[i, ], as.matrix(covEst[i, , ]),
                      nuEst[i], log = TRUE)
    }
    lli    <- apply(lli, 2, pmax, -30)
    lnresp <- matrixStats::rowLogSumExps(lli)
    lliCur <- sum(lnresp)
    traceLli <- c(traceLli, lliCur)
    if (iter > 1)
      Dlli <- (traceLli[iter + 1] - traceLli[iter - 1]) / 2
    if (Dlli < tol) break

    lresp <- sweep(lli, 1, lnresp, "-")
    resp  <- exp(lresp)

    for (i in seq(Nc)) {
      Xc_i <- sweep(X, 2L, muEst[i, ])
      Sigma_i <- as.matrix(covEst[i, , ])
      maha <- tryCatch(mahalanobis(Xc_i, FALSE, Sigma_i),
                       error = function(e) {
                         L <- t(chol(Sigma_i +
                                     diag(5 * .Machine$double.eps, Nd, Nd)))
                         colSums(forwardsolve(L, t(Xc_i))^2)
                       })
      lU[, i] <- log(nuEst[i] + Nd) - log(nuEst[i] + maha)
    }

    lzSum <- colLogSumExps(lresp)
    zSum  <- exp(lzSum)
    pEst  <- exp(lzSum - logN); pEst <- pEst / sum(pEst)

    lW    <- lresp + lU
    wSum  <- exp(colLogSumExps(lW))
    for (i in seq(Nc)) {
      xc <- exp(0.5 * (lW[, i] - lzSum[i])) *
            sweep(X, 2, muEst[i, ], check.margin = FALSE)
      covEst[i, , ] <- crossprod(xc)
      muEst[i, ]    <- colSums(exp(lW[, i]) * X) / wSum[i]
      for (j in 1:Nd)
        covEst[i, j, j] <- max(covEst[i, j, j], .Machine$double.eps)
    }

    c1 <- 1 + colSums(resp * (lU - exp(lU))) / zSum
    for (i in seq(Nc))
      nuEst[i] <- optimize(nu_ml(c1[i]), c(0, 150))$minimum

    iter <- iter + 1
  }

  o      <- order(pEst, decreasing = TRUE)
  pEst   <- pEst[o]
  muEst  <- muEst[o, , drop = FALSE]
  covEst <- covEst[o, , , drop = FALSE]
  nuEst  <- nuEst[o, drop = FALSE]

  invisible(list(p = pEst, center = muEst, cov = covEst, nu = nuEst))
}


## ------------------------------------------------------------
## EM_mnmm -- multivariate normal mixture EM (Student-t init).
## Returns an mvnormMix object with EM tracing attributes
## (lli, df, nobs, traceLli, traceMix, x) used by logLik.EM.
## ------------------------------------------------------------

## Helper: pack (mean, cov) of one component into a vector of
## (mean, sd, lower-tri correlations).  1-D fast path is just (m, s).
mv2vec <- function(mean, sigma) {
  Nd <- length(mean)
  if (Nd != 1) c(mean, sqrt(diag(sigma)), cov2cor(sigma)[lower.tri(sigma)])
  else         c(mean, sqrt(sigma))
}

EM_mnmm <- function(X, Nc, mix_init, Ninit = 50, verbose = FALSE,
                    Niter.max = 500, tol, Neps,
                    eps = c(w = 0.005, m = 0.005, s = 0.005)) {

  if (!is.matrix(X)) X <- matrix(X, ncol = 1)
  N  <- dim(X)[1]
  Nd <- dim(X)[2]
  assert_that(N + Nc >= Ninit)

  if (missing(mix_init)) {
    mix_init <- EM_msmm(X, Nc, Ninit = Ninit, verbose = verbose,
                        Niter.max = round(Niter.max / 2), tol = 0.1)
  }
  pEst   <- mix_init$p
  muEst  <- mix_init$center
  covEst <- mix_init$cov

  est2par <- function(p, mu, cov) {
    est <- rbind(logit(p),
                 matrix(sapply(1:Nc,
                               function(i) mv2vec(mu[i, ], cov[i, , ])),
                        ncol = Nc))
    est[(1 + Nd + 1):(1 + 2 * Nd), ] <- log(est[(1 + Nd + 1):(1 + 2 * Nd), ])
    est
  }

  if (missing(tol))  { checkTol <- FALSE; tol <- -1 } else checkTol <- TRUE
  if (missing(Neps)) { checkEps <- FALSE; Neps <- 5 } else checkEps <- TRUE
  if (!checkTol & !checkEps) checkEps <- TRUE
  assert_that(Neps > 1, ceiling(Neps) == floor(Neps))

  if (length(eps) == 1) eps <- rep(10^(-eps), 3)

  cov.df  <- (Nd - 1) * Nd / 2 + Nd
  df      <- Nc * (Nd + cov.df) + Nc - 1
  df.comp <- cov.df + Nd + 1

  eps <- c(eps[1], rep(eps[2], Nd), rep(eps[3], Nd),
           rep(eps[2], (Nd - 1) * Nd / 2))

  iter      <- 0
  logN      <- log(N)
  traceMix  <- list()
  traceLli  <- c()
  Dlli      <- Inf
  runMixPar <- array(-Inf, dim = c(Neps, df.comp, Nc))
  runOrder  <- 0:(Neps - 1)
  Npar      <- if (Nc == 1) df else df + 1
  lli       <- array(-20, dim = c(N, Nc))

  while (iter < Niter.max) {
    for (i in seq(Nc)) {
      lli[, i] <- log(pEst[i]) +
        dmvnorm(X, muEst[i, ], as.matrix(covEst[i, , ]),
                log = TRUE, checkSymmetry = FALSE)
    }
    lli    <- apply(lli, 2, pmax, -30)
    lnresp <- matrixStats::rowLogSumExps(lli)
    lliCur <- sum(lnresp)
    traceMix <- c(traceMix, list(list(p = pEst, mean = muEst, sigma = covEst)))
    traceLli <- c(traceLli, lliCur)
    if (iter > 1)
      Dlli <- (traceLli[iter + 1] - traceLli[iter - 1]) / 2

    if (Nc > 1) {
      smean <- apply(runMixPar[order(runOrder), , , drop = FALSE], c(2, 3),
                     function(x) mean(abs(diff(x))))
      eps.converged <- sum(sweep(smean, 1, eps, "-") < 0)
    } else {
      smean <- apply(runMixPar[order(runOrder), -1, , drop = FALSE], c(2, 3),
                     function(x) mean(abs(diff(x))))
      eps.converged <- sum(sweep(smean, 1, eps[-1], "-") < 0)
    }
    if (is.na(eps.converged)) eps.converged <- 0
    if (checkTol & Dlli < tol) break
    if (iter >= Neps & checkEps & eps.converged == Npar) break

    lresp <- sweep(lli, 1, lnresp, "-")
    lzSum <- colLogSumExps(lresp)
    pEst  <- exp(lzSum - logN); pEst <- pEst / sum(pEst)

    for (i in seq(Nc)) {
      upd <- cov.wt(X, exp(lresp[, i] - lzSum[i]), method = "ML")
      muEst[i, ]    <- upd$center
      covEst[i, , ] <- upd$cov
      for (j in 1:Nd)
        covEst[i, j, j] <- max(covEst[i, j, j], .Machine$double.eps)
    }

    ind <- 1 + iter %% Neps
    runMixPar[ind, , ] <- est2par(pEst, muEst, covEst)
    runOrder[ind] <- iter
    iter <- iter + 1
  }
  if (iter + 1 == Niter.max) warning("Maximum number of iterations reached.")

  o      <- order(pEst, decreasing = TRUE)
  pEst   <- pEst[o]
  muEst  <- muEst[o, , drop = FALSE]
  covEst <- covEst[o, , , drop = FALSE]

  mixEst <- do.call(mixmvnorm,
    lapply(1:Nc, function(i)
      c(pEst[i], muEst[i, , drop = TRUE], matrix(covEst[i, , ], Nd, Nd))))
  attr(mixEst, "df")       <- df
  attr(mixEst, "nobs")     <- N
  attr(mixEst, "lli")      <- lliCur
  attr(mixEst, "Nc")       <- Nc
  attr(mixEst, "tol")      <- tol
  attr(mixEst, "traceLli") <- traceLli
  attr(mixEst, "traceMix") <- lapply(traceMix, function(est)
    suppressWarnings(do.call(mixmvnorm,
      lapply(1:Nc, function(i)
        c(est$p[i], est$mean[i, , drop = FALSE],
          matrix(est$sigma[i, , ], Nd, Nd))))))
  attr(mixEst, "x")        <- X
  class(mixEst) <- c("EM", "EMmvnmm", "mvnormMix", "mix")
  mixEst
}


## ------------------------------------------------------------
## EM_nmm -- univariate wrapper called by mixfit(type = "norm")
## ------------------------------------------------------------

EM_nmm <- function(X, Nc, mix_init, verbose = FALSE, Niter.max = 500,
                   tol, Neps, eps = c(w = 0.005, m = 0.005, s = 0.005)) {
  if (is.matrix(X)) assert_matrix(X, any.missing = FALSE, ncols = 1)
  mixEst <- EM_mnmm(X = X, Nc = Nc, mix_init = mix_init,
                    verbose = verbose, Niter.max = Niter.max,
                    tol = tol, Neps = Neps, eps = eps)
  rownames(mixEst) <- c("w", "m", "s")
  class(mixEst) <- c("EM", "EMnmm", "normMix", "mix")
  attr(mixEst, "traceMix") <- lapply(attr(mixEst, "traceMix"), function(x) {
    class(x) <- class(mixEst); rownames(x) <- rownames(mixEst); x
  })
  likelihood(mixEst) <- "normal"
  mixEst
}
