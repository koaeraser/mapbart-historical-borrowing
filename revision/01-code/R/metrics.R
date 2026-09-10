# metrics.R: per-replicate metrics and their aggregation over replicates.
#
# Metric definitions follow Table "ate_metrics" of paper/tables_figures.tex.
# For an estimand with posterior draws tau_b (b = 1..B), truth tau_true and
# null value tau_0 (0 for differences, 1 for ratios):
#   pm        posterior mean
#   bias      pm - tau_true                                (averaged over r)
#   sd        posterior SD                                  (averaged over r)
#   rmse_draws sqrt(B^-1 sum_b (tau_b - tau_true)^2), the table's formula
#             (averaged over r)
#   rmse      sqrt(mean_r (pm_r - tau_true)^2), the Monte Carlo RMSE of the
#             posterior mean. The paper's tables are consistent with this
#             definition and not with rmse_draws (e.g. BART-PP Sc1 reports
#             SD 0.21 and RMSE 0.20, which the draw-based formula cannot
#             produce since it is bounded below by the posterior SD). Both
#             are reported; "rmse" is the column matched to the paper.
#   cover     1{tau_true in the central 95% credible interval}
#   reject    1{tau_0 not in the central 95% credible interval} ("power" under
#             the alternative, Type I error rate under the null; the paper's
#             stated definition)
#   reject90  the same with the central 90% interval (the definition the
#             published power column matches; reports/reference_check.md)
# Arm-specific metrics use the same definitions for the trt and ctrl draws.
#
# Additional metrics for the revision (NA when the method does not provide
# the inputs):
#   surf_rmse, surf_rmse_R, surf_rmse_Rc  RMSE of the posterior-mean control
#             surface against the true RCT control surface at the RCT
#             profiles, overall and restricted to profiles in R and in R^c
#   map_R, map_Rc, map_all  mean of the method's P(borrow) at the RCT control
#             profiles in R, in R^c, and overall
#   ess_prior, ess_realized  as returned by the method.
#   mapspike_R, mapspike_Rc, g_mean_R, g_mean_Rc, ess_ceiling, ess_capped
#             secondary LRC-BART diagnostics read from fit$extra / fit$ess
#             (all-spike map, mean posterior g by region, ESS ceiling, cap
#             flag; ess_capped aggregates to the fraction of capped replicates).

suppressPackageStartupMessages(library(data.table))

estimand_groups <- function(outcome) {
  if (outcome == "gaussian") {
    list(list(ate = "ate", trt = "trt", ctrl = "ctrl", suffix = "", null = 0))
  } else {
    list(list(ate = "ate_med", trt = "trt_med", ctrl = "ctrl_med", suffix = "_med", null = 1),
         list(ate = "ate_rmst", trt = "trt_rmst", ctrl = "ctrl_rmst", suffix = "_rmst", null = 1))
  }
}

one_estimand_metrics <- function(draws, truth, null) {
  draws <- draws[is.finite(draws)]
  q <- quantile(draws, c(0.025, 0.975), names = FALSE)
  q90 <- quantile(draws, c(0.05, 0.95), names = FALSE)
  pm <- mean(draws)
  c(pm = pm, truth = truth, bias = pm - truth, sd = sd(draws),
    rmse_draws = sqrt(mean((draws - truth)^2)),
    lo = q[1], hi = q[2],
    cover = as.numeric(truth >= q[1] & truth <= q[2]),
    reject = as.numeric(null < q[1] | null > q[2]),
    reject90 = as.numeric(null < q90[1] | null > q90[2]))
}

# fit: the object returned by a method; dat: the gen_data() output.
compute_metrics <- function(fit, dat) {
  outcome <- dat$constants$outcome
  out <- list()
  for (g in estimand_groups(outcome)) {
    for (arm in c("ate", "trt", "ctrl")) {
      nm <- g[[arm]]
      m <- one_estimand_metrics(fit$draws[[nm]], dat$truth[[nm]],
                                null = if (arm == "ate") g$null else NA)
      names(m) <- paste0(arm, "_", names(m), g$suffix)
      out[[nm]] <- m
    }
  }
  res <- as.list(unlist(unname(out)))
  # Control-surface RMSE at the RCT profiles (overall, in R, in R^c).
  if (!is.null(fit$surface) && length(fit$surface) == length(dat$mu0_rct)) {
    err <- fit$surface - dat$mu0_rct
    R <- dat$region_rct
    res$surf_rmse <- sqrt(mean(err^2))
    res$surf_rmse_R <- if (any(R)) sqrt(mean(err[R]^2)) else NA_real_
    res$surf_rmse_Rc <- if (any(!R)) sqrt(mean(err[!R]^2)) else NA_real_
  } else {
    res$surf_rmse <- res$surf_rmse_R <- res$surf_rmse_Rc <- NA_real_
  }
  # Borrowing map summaries at the RCT control profiles.
  Rc <- dat$region_rct[dat$z_rct == 0]
  if (!is.null(fit$map) && length(fit$map) == length(Rc)) {
    res$map_all <- mean(fit$map)
    res$map_R <- if (any(Rc)) mean(fit$map[Rc]) else NA_real_
    res$map_Rc <- if (any(!Rc)) mean(fit$map[!Rc]) else NA_real_
  } else {
    res$map_all <- res$map_R <- res$map_Rc <- NA_real_
  }
  res$ess_prior <- if (!is.null(fit$ess$prior)) as.numeric(fit$ess$prior) else NA_real_
  res$ess_realized <- if (!is.null(fit$ess$realized)) as.numeric(fit$ess$realized) else NA_real_
  # Secondary map diagnostics carried in fit$extra (LRC-BART): the ensemble-level
  # all-spike map, the mean posterior discrepancy g(x) by region, the ESS
  # ceiling and the feasibility cap flag.
  region_mean <- function(v, sel) if (!is.null(v) && length(v) == length(Rc) && any(sel)) mean(v[sel]) else NA_real_
  res$mapspike_R <- region_mean(fit$extra$map_spike, Rc)
  res$mapspike_Rc <- region_mean(fit$extra$map_spike, !Rc)
  res$g_mean_R <- region_mean(fit$extra$g_mean, Rc)
  res$g_mean_Rc <- region_mean(fit$extra$g_mean, !Rc)
  res$ess_ceiling <- if (!is.null(fit$ess$ceiling)) as.numeric(fit$ess$ceiling) else NA_real_
  res$ess_capped <- if (!is.null(fit$extra$ess_capped)) as.numeric(fit$extra$ess_capped) else NA_real_
  as.data.frame(res)
}

# Aggregate a data.table of per-replicate rows (with columns scenario, method,
# rep) into the Table 1 layout. Returns a data.table with one row per
# scenario x method.
aggregate_metrics <- function(dt, outcome) {
  dt <- as.data.table(dt)
  groups <- estimand_groups(outcome)
  agg_one <- function(d) {
    res <- list(n_rep = nrow(d))
    for (g in groups) {
      sfx <- g$suffix
      for (arm in c("ate", "trt", "ctrl")) {
        pm <- d[[paste0(arm, "_pm", sfx)]]; tr <- d[[paste0(arm, "_truth", sfx)]]
        res[[paste0(arm, "_bias", sfx)]] <- mean(pm - tr)
        res[[paste0(arm, "_sd", sfx)]] <- mean(d[[paste0(arm, "_sd", sfx)]])
        res[[paste0(arm, "_rmse", sfx)]] <- sqrt(mean((pm - tr)^2))
        res[[paste0(arm, "_rmse_draws", sfx)]] <- mean(d[[paste0(arm, "_rmse_draws", sfx)]])
        res[[paste0(arm, "_cover", sfx)]] <- mean(d[[paste0(arm, "_cover", sfx)]])
        if (arm == "ate") {
          res[[paste0("ate_power", sfx)]] <- mean(d[[paste0("ate_reject", sfx)]])
          r90 <- paste0("ate_reject90", sfx)
          res[[paste0("ate_power90", sfx)]] <- if (r90 %in% names(d)) mean(d[[r90]]) else NA_real_
        }
      }
    }
    for (v in c("surf_rmse", "surf_rmse_R", "surf_rmse_Rc", "map_all", "map_R", "map_Rc",
                "ess_prior", "ess_realized", "mapspike_R", "mapspike_Rc", "g_mean_R", "g_mean_Rc",
                "ess_ceiling", "ess_capped", "seconds")) {
      res[[v]] <- if (v %in% names(d)) mean(d[[v]], na.rm = TRUE) else NA_real_
    }
    as.data.table(res)
  }
  dt[, agg_one(.SD), by = .(scenario, method)]
}

# Compact Table 1 view: scenario, method, ATE bias/SD/RMSE/coverage/power and
# the arm-specific bias/RMSE, for the primary estimand.
table1_view <- function(agg, outcome) {
  sfx <- if (outcome == "gaussian") "" else "_med"
  cols <- c("scenario", "method", "n_rep",
            paste0("ate_", c("bias", "sd", "rmse", "cover", "power", "power90"), sfx),
            paste0("trt_", c("bias", "rmse"), sfx),
            paste0("ctrl_", c("bias", "rmse"), sfx))
  out <- as.data.frame(agg)[, cols]
  names(out) <- c("scenario", "method", "n_rep", "bias", "sd", "rmse", "cov", "pow", "pow90",
                  "trt_bias", "trt_rmse", "ctrl_bias", "ctrl_rmse")
  out
}
