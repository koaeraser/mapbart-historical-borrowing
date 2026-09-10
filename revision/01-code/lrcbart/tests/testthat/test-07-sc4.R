# Test 7 (spec test 3): H_f = 50, H_g at the package default (5 after Phase 2b) on a leaf-heterogeneous DGP with a
# delta = 2 regional shift in the RWD. The borrowing map must separate the
# regions and the standardized ATE must be close to the truth.
sc4_data <- function(seed, n_rct = 300, n_rwd = 300, delta = 2, sig1 = 1.0, sig2 = 1.2) {
  set.seed(seed)
  mu <- c(1, -1, 0.5, -1, 2, -2, 2, -3, 1, 0)
  gen_x <- function(n) sweep(matrix(rnorm(n * 10), n, 10), 2, mu, "+")
  fx <- function(X) -0.5 * X[, 1]^2 + 0.5 * X[, 2] - 0.5 * abs(X[, 3] - 0.5) - 1.35 * tanh(X[, 5] - 2) - 0.8 * tanh(X[, 6] + 2)
  x_rct <- gen_x(n_rct); x_rwd <- gen_x(n_rwd)
  trt <- rbinom(n_rct, 1, 2 / 3)
  y_rct <- fx(x_rct) + 1 * trt + rnorm(n_rct, 0, sig1)
  inR_rwd <- x_rwd[, 5] > 2; inR_rct <- x_rct[, 5] > 2
  y_rwd <- fx(x_rwd) + delta * (!inR_rwd) + rnorm(n_rwd, 0, sig2)
  list(x_trt = x_rct[trt == 1, ], y_trt = y_rct[trt == 1],
       x_ctrl = x_rct[trt == 0, ], y_ctrl = y_rct[trt == 0],
       x_rwd = x_rwd, y_rwd = y_rwd, x_rct = x_rct, inR_ctrl = inR_rct[trt == 0],
       ate_true = 1, mu0_ctrl = fx(x_rct[trt == 0, ]))
}

test_that("Sc4-type DGP: borrowing map separates regions, ATE bias small", {
  seeds <- c(101, 102, 103)
  res <- t(sapply(seeds, function(s) {
    d <- sc4_data(s)
    cal <- lrc_ess_calibrate(d$y_rwd, d$x_rwd, d$x_rct, N_target = 100, n_burn = 500, n_draw = 1000, seed = s,
                             y_rct_ctrl = d$y_ctrl, x_rct_ctrl = d$x_ctrl)
    fit <- lrc_bart(c(d$y_ctrl, d$y_rwd), rbind(d$x_ctrl, d$x_rwd),
                    c(rep(1, length(d$y_ctrl)), rep(2, length(d$y_rwd))),
                    s0_sq = cal$s0_sq, calib = cal, n_burn = 1000, n_draw = 1000, seed = s + 1)
    ftrt <- bart_fit(d$y_trt, d$x_trt, n_burn = 1000, n_draw = 1000, seed = s + 2)
    ate <- lrc_ate(ftrt, fit, d$x_rct, "gaussian")
    bm <- lrc_borrowing_map(fit, d$x_ctrl, type = "abs_g", c = 0.5)$mean
    bs <- lrc_borrowing_map(fit, d$x_ctrl, type = "all_spike")$mean
    gm <- colMeans(lrc_discrepancy(fit, d$x_ctrl))
    mu0 <- colMeans(predict(fit, d$x_ctrl))
    c(s0 = cal$s0_sq, ess0 = cal$ess0_at_s0, ceiling = cal$ceiling, capped = cal$ess_capped,
      ess_real = lrc_ess_realized(fit, d$x_ctrl)$mean,
      B_R = mean(bm[d$inR_ctrl]), B_Rc = mean(bm[!d$inR_ctrl]),
      Bspike_R = mean(bs[d$inR_ctrl]), Bspike_Rc = mean(bs[!d$inR_ctrl]),
      g_R = mean(gm[d$inR_ctrl]), g_Rc = mean(gm[!d$inR_ctrl]),
      bias_R = mean(mu0[d$inR_ctrl] - d$mu0_ctrl[d$inR_ctrl]),
      bias_Rc = mean(mu0[!d$inR_ctrl] - d$mu0_ctrl[!d$inR_ctrl]),
      ate = mean(ate$ate), ate_sd = sd(ate$ate), spike_frac = mean(fit$K0 / fit$L_g))
  }))
  rownames(res) <- paste0("seed", seeds)
  cat("\n[T7] per-replicate summaries (truth: ATE = 1, g = 0 in R, g = -2 in R^c; B = P(|g| < 0.5), Bspike = all-spike map):\n")
  print(round(res, 3))
  cat(sprintf("[T7] mean B_R = %.3f (>= 0.7), mean B_Rc = %.3f (<= 0.15), mean ATE bias = %.3f (< 0.1), mean |ctrl bias| R^c = %.3f\n",
              mean(res[, "B_R"]), mean(res[, "B_Rc"]), mean(res[, "ate"]) - 1, mean(abs(res[, "bias_Rc"]))))
  expect_lte(mean(res[, "B_Rc"]), 0.15)
  expect_lt(abs(mean(res[, "ate"]) - 1), 0.1)
  expect_lt(mean(abs(res[, "bias_Rc"])), 0.2)
  expect_lt(abs(mean(res[, "g_Rc"]) + 2), 0.35)
  expect_lt(abs(mean(res[, "g_R"])), 0.25)
  # compatibility map P(|g(x)| < 0.5 | data): >= 0.7 in R (post-verification
  # criterion; the all-spike map is a diagnostic only, see NOTES.md)
  expect_gte(mean(res[, "B_R"]), 0.7)
})
