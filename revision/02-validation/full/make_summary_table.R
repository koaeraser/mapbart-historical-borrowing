# make_summary_table(): the paper's table layout from summarize_sim()'s full
# aggregate. One row per scenario x method (x estimand for survival):
# bias, sd, rmse, cov, pow (95% interval excludes the null, the paper's stated
# definition), pow90 (90% interval), arm-specific bias and RMSE, the Sc4
# control-surface RMSE overall / in R / in R^c, the compatibility-map means
# (map_R, map_Rc, map_all), the mean posterior g by region, and the ESS
# columns (prior ESS_0, realized, ceiling, capped fraction). Also n_rep and
# mean seconds per fit.
make_summary_table <- function(agg, outcome, extra_cols = NULL) {
  agg <- as.data.frame(agg)
  sfxs <- if (outcome == "gaussian") c("") else c("_med", "_rmst")
  out <- list()
  for (sfx in sfxs) {
    cols <- c(paste0("ate_", c("bias", "sd", "rmse", "cover", "power", "power90"), sfx),
              paste0("trt_", c("bias", "rmse"), sfx), paste0("ctrl_", c("bias", "rmse"), sfx))
    d <- agg[, c("scenario", "method", "n_rep", cols)]
    names(d) <- c("scenario", "method", "n_rep", "bias", "sd", "rmse", "cov", "pow", "pow90",
                  "trt_bias", "trt_rmse", "ctrl_bias", "ctrl_rmse")
    d$estimand <- if (sfx == "") "ate" else sub("_", "", sfx)
    diag <- c("surf_rmse", "surf_rmse_R", "surf_rmse_Rc", "map_R", "map_Rc", "map_all",
              "mapspike_R", "mapspike_Rc", "g_mean_R", "g_mean_Rc",
              "ess_prior", "ess_realized", "ess_ceiling", "ess_capped", "seconds")
    d <- cbind(d, agg[, intersect(diag, names(agg))])
    out[[length(out) + 1]] <- d
  }
  out <- do.call(rbind, out)
  if (!is.null(extra_cols)) out <- cbind(extra_cols, out)
  out[order(out$estimand, out$scenario, out$method), ]
}
