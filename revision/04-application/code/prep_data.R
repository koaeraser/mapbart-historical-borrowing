# prep_data.R: load the de-identified merged EloKRd + UCMM analysis data
# (yunxuan-repo/mapbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n283.RData)
# and build the covariate matrices. Two codings are returned:
#   X_orig: Yunxuan's prep_covariates() (K-1 dummies per factor, levels sorted
#           alphabetically), reproduced verbatim so that the paper's reference
#           rows can be rerun. Because ethnicity is labelled "Non-Hispanic" in
#           EloKRd and "Not Hispanic or Latino" in UCMM, this coding gives two
#           non-Hispanic dummies, one per source (see data_check.md).
#   X_harm: harmonized coding, 7 columns: age, male, race_Black,
#           race_Other (Asian and Other/Unknown pooled, as in the paper's Table
#           baseline), hispanic, high_risk_cyto, asct.
# Times are converted from months to years, as in every script of the repo.
ROOT <- normalizePath(file.path(Sys.getenv("LRC_REVISION_ROOT", "revision"), ".."))
DATA_FILE <- file.path(ROOT, "mapbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n283.RData")
COV_COLS <- c("age", "sex", "race", "ethnicity", "high_risk_cyto", "transplant_off_protocol")

prep_covariates_orig <- function(df, cols) {
  X <- data.frame(row.names = 1:nrow(df))
  for (v in cols) {
    vals <- df[[v]]
    if (is.numeric(vals)) {
      X[[v]] <- vals
    } else {
      vals <- as.factor(vals); lvls <- levels(vals)
      if (length(lvls) == 2) X[[v]] <- as.integer(vals == lvls[2])
      else for (k in 2:length(lvls))
        X[[paste0(v, "_", gsub("[/ ]", "_", lvls[k]))]] <- as.integer(vals == lvls[k])
    }
  }
  X
}

prep_covariates_harm <- function(df) {
  race3 <- ifelse(df$race == "White", "White", ifelse(df$race == "Black", "Black", "Other"))
  data.frame(age = df$age,
             male = as.integer(df$sex == "Male"),
             race_Black = as.integer(race3 == "Black"),
             race_Other = as.integer(race3 == "Other"),
             hispanic = as.integer(df$ethnicity == "Hispanic or Latino"),
             high_risk_cyto = as.integer(df$high_risk_cyto),
             asct = as.integer(df$transplant_off_protocol))
}

load_mm <- function(outcome = c("PFS", "OS")) {
  outcome <- match.arg(outcome)
  e <- new.env(); load(DATA_FILE, envir = e); d <- e$merged
  stopifnot(nrow(d) == 283, sum(d$trt == 1) == 30, sum(d$trt == 0) == 253)
  time <- if (outcome == "PFS") d$pfs_months / 12 else d$os_months / 12
  status <- if (outcome == "PFS") d$pfs_status else d$os_status
  list(outcome = outcome, time = time, status = as.integer(status), trt = d$trt,
       trt_idx = which(d$trt == 1), ctrl_idx = which(d$trt == 0),
       X_orig = prep_covariates_orig(d, COV_COLS),
       X_harm = prep_covariates_harm(d))
}
