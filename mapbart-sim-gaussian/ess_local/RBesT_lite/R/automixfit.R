## Trimmed: gMAPpred branch removed (we always pass a gMAP object);
## verbose branch consolidated.

automixfit <- function(sample, Nc = seq(1, 4), k = 6, thresh = -Inf,
                       verbose = FALSE, ...) {
  assert_that(all(diff(Nc) >= 1))
  models <- list(); aic <- c(); ic <- Inf
  for (i in seq_along(Nc)) {
    icLast <- ic
    run <- if (verbose) mixfit(sample, Nc = Nc[i], verbose = verbose, ...)
           else suppressMessages(mixfit(sample, Nc = Nc[i],
                                        verbose = verbose, ...))
    ic     <- AIC(run, k = k)
    aic    <- c(aic, ic)
    models <- c(models, list(run))
    if (icLast - ic < thresh) break
  }
  names(models) <- Nc[1:length(models)]
  models <- models[order(aic)]
  bestfit <- models[[1]]
  attr(bestfit, "models") <- models
  bestfit
}
