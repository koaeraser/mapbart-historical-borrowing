## Load a self-contained "lite" copy of RBesT into the current
## environment.  No need to install or library(RBesT) -- the R sources
## live alongside this file in ./R, and the Stan model lives in ./stan
## and is compiled on first use (cached by rstan thereafter).
##
## After sourcing, the symbols gMAP, automixfit, ess (and friends) are
## defined in .GlobalEnv exactly as the package would expose them.

local({
  ## --- CRAN deps actually used along the gMAP -> automixfit -> ess
  ## pipeline.  We attach (not just loadNamespace) the ones whose S3
  ## generics or operators the package code calls without the pkg::
  ## prefix.
  needed <- c(
    "methods", "stats", "utils",
    "Rcpp", "RcppParallel", "rstan", "rstantools",
    "assertthat", "checkmate", "Formula",
    "matrixStats", "abind", "mvtnorm",
    "dplyr", "ggplot2", "rlang", "jsonlite", "lifecycle",
    "bayesplot", "posterior"
  )
  for (p in needed) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, quietly = TRUE)
    )
  }
})

## Resolve our own directory (works under Rscript and source()).
## Allow override by pre-defining `.RBesT_lite_dir` or setting env var
## RBESLITE_DIR, useful when this file is sourced from another script.
.RBesT_lite_dir <- local({
  override <- Sys.getenv("RBESLITE_DIR", unset = NA)
  if (!is.na(override) && nzchar(override) && dir.exists(override))
    return(normalizePath(override))
  if (exists(".RBesT_lite_dir", envir = .GlobalEnv) &&
      dir.exists(.GlobalEnv$.RBesT_lite_dir))
    return(normalizePath(.GlobalEnv$.RBesT_lite_dir))
  this <- tryCatch(
    sys.frames()[[1]]$ofile, error = function(e) NULL
  )
  if (is.null(this) || !nzchar(this)) {
    for (fr in rev(sys.frames())) {
      ofile <- fr$ofile
      if (!is.null(ofile) && nzchar(ofile) &&
          basename(ofile) == "load.R") { this <- ofile; break }
    }
  }
  if (is.null(this) || !nzchar(this)) {
    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    if (length(fa)) this <- sub("^--file=", "", fa[1])
  }
  if (is.null(this) || !nzchar(this))
    this <- normalizePath("RBesT_lite/load.R", mustWork = FALSE)
  dirname(normalizePath(this, mustWork = FALSE))
})

## --- internal package data (calibration_data, calibration_meta,
## pkg_create_date, pkg_sha) used by a couple of helper functions.
load(file.path(.RBesT_lite_dir, "R", "sysdata.rda"), envir = .GlobalEnv)

## --- source every R file we kept from the package.  Order doesn't
## matter for `<-` assignments at the top level.
.RBesT_lite_files <- list.files(file.path(.RBesT_lite_dir, "R"),
                                pattern = "\\.R$", full.names = TRUE)
for (.f in .RBesT_lite_files) {
  source(.f, local = FALSE, keep.source = FALSE)
}
rm(.f, .RBesT_lite_files)

## --- Stan models.  Three variants of the MAP model:
##   standard  -- theta_h = mu + tau xi_h, mu ~ N(0, sd^2)  (default MAP)
##   hybrid    -- theta_h = tau xi_h        (mu fixed at 0)  (two-arm hybrid)
##   synthetic -- theta_h ~ N(0, sd^2), theta_star ~ N(theta_h, tau^2)
##                                              (single-arm synthetic control)
## RBesT normally ships precompiled C++ modules; we compile from
## source on first use and rely on rstan's .rds cache thereafter.
rstan::rstan_options(auto_write = TRUE)
.RBesT_lite_versions <- c("standard", "hybrid", "synthetic")
if (!exists("stanmodels", envir = .GlobalEnv) ||
    !is.list(.GlobalEnv$stanmodels) ||
    !all(.RBesT_lite_versions %in% names(.GlobalEnv$stanmodels))) {
  .stanmodels <- lapply(.RBesT_lite_versions, function(v) {
    f <- file.path(.RBesT_lite_dir, "stan", paste0("gMAP_", v, ".stan"))
    message("RBesT_lite: compiling gMAP_", v,
            ".stan (one-time, cached by rstan) ...")
    rstan::stan_model(file = f, model_name = paste0("gMAP_", v))
  })
  names(.stanmodels) <- .RBesT_lite_versions
  assign("stanmodels", .stanmodels, envir = .GlobalEnv)
  rm(.stanmodels)
}
rm(.RBesT_lite_versions)

invisible(NULL)
