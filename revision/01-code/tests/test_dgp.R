# test_dgp.R: checks that the data-generating mechanism reproduces the
# paper's summary facts. Run from the 01-code root:
#   Rscript tests/test_dgp.R
suppressPackageStartupMessages(library(testthat))
CODE_ROOT <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
source(file.path(CODE_ROOT, "R", "dgp.R"))

n_rep <- 50
seeds <- seq_len(n_rep)

test_that("Sc1 to Sc3 Gaussian: sizes, allocation, true ATE", {
  for (sc in c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5")) {
    n_ctrl <- ate <- ate_pop <- numeric(n_rep)
    for (r in seeds) {
      d <- gen_data_id(sc, "gaussian", seed = r, n_mc = 2e4)
      expect_equal(nrow(d$rct_trt) + nrow(d$rct_ctrl), 300)
      expect_equal(nrow(d$rwd_ctrl), 300)
      n_ctrl[r] <- nrow(d$rct_ctrl)
      ate[r] <- d$truth$ate
      ate_pop[r] <- d$truth_pop$ate
    }
    # 2:1 allocation: control share around 100 of 300 (Binomial(300, 1/3))
    expect_true(abs(mean(n_ctrl) - 100) < 3, label = paste(sc, "control share"))
    # true Gaussian ATE = delta = 1 exactly (no effect modification), both
    # sample-standardized and population versions
    expect_true(all(abs(ate - 1) < 1e-10), label = paste(sc, "true ATE"))
    expect_true(all(abs(ate_pop - 1) < 1e-10), label = paste(sc, "true ATE pop"))
  }
})

test_that("Sc2 and Sc3: source imbalance in the confounders", {
  diff_sc1 <- diff_sc2 <- diff_u <- diff_x_rho0 <- diff_x_rho5 <- numeric(n_rep)
  for (r in seeds) {
    d1 <- gen_data_id("Sc1", "gaussian", seed = r, n_mc = 0)
    d2 <- gen_data_id("Sc2", "gaussian", seed = r, n_mc = 0)
    d3 <- gen_data_id("Sc3_rho0", "gaussian", seed = r, n_mc = 0)
    d5 <- gen_data_id("Sc3_rho0.5", "gaussian", seed = r, n_mc = 0)
    diff_sc1[r] <- mean(d1$x_rct[, 5] + d1$x_rct[, 6]) - mean(d1$rwd_ctrl$X5 + d1$rwd_ctrl$X6)
    diff_sc2[r] <- mean(d2$x_rct[, 5] + d2$x_rct[, 6]) - mean(d2$rwd_ctrl$X5 + d2$rwd_ctrl$X6)
    diff_u[r] <- mean(rowSums(d3$latent$U_rct)) - mean(rowSums(d3$latent$U_rwd))
    diff_x_rho0[r] <- mean(d3$x_rct[, 5] + d3$x_rct[, 6]) - mean(d3$rwd_ctrl$X5 + d3$rwd_ctrl$X6)
    diff_x_rho5[r] <- mean(d5$x_rct[, 5] + d5$x_rct[, 6]) - mean(d5$rwd_ctrl$X5 + d5$rwd_ctrl$X6)
  }
  # Sc1: no imbalance; Sc2: RCT has lower X5 + X6 (beta_D = -1.2); Sc3: the
  # imbalance is in U (rho = 0 leaves X balanced, rho = 0.5 makes X a proxy)
  expect_true(abs(mean(diff_sc1)) < 0.05)
  expect_true(mean(diff_sc2) < -0.8)
  expect_true(mean(diff_u) < -0.8)
  expect_true(abs(mean(diff_x_rho0)) < 0.05)
  expect_true(mean(diff_x_rho5) < -0.3)
})

test_that("Sc1 to Sc3 survival: 3:1 allocation, median ratio 1.4, censoring", {
  for (sc in c("Sc1", "Sc2", "Sc3_rho0")) {
    n_ctrl <- ev_rct <- ev_rwd <- mr <- mr_pop <- rr <- numeric(n_rep)
    for (r in seeds) {
      d <- gen_data_id(sc, "survival", seed = r, n_mc = 2e4)
      n_ctrl[r] <- nrow(d$rct_ctrl)
      ev_rct[r] <- mean(c(d$rct_trt$status, d$rct_ctrl$status))
      ev_rwd[r] <- mean(d$rwd_ctrl$status)
      mr[r] <- d$truth$ate_med; mr_pop[r] <- d$truth_pop$ate_med
      rr[r] <- d$truth$ate_rmst
      expect_true(all(d$rct_trt$time > 0 & d$rct_trt$time <= 3))
    }
    expect_true(abs(mean(n_ctrl) - 75) < 3, label = paste(sc, "3:1 allocation"))
    expect_true(all(abs(mr - 1.4) < 1e-6), label = paste(sc, "median ratio"))
    expect_true(all(abs(mr_pop - 1.4) < 1e-6), label = paste(sc, "median ratio pop"))
    expect_true(all(rr > 1 & rr < 1.4), label = paste(sc, "RMST ratio"))
    # event proportions in a plausible range (dropout Exp(0.02) is light;
    # administrative censoring at 3 - U(0,1) dominates)
    expect_true(mean(ev_rct) > 0.5 & mean(ev_rct) < 0.95, label = paste(sc, "RCT events"))
    expect_true(mean(ev_rwd) > 0.5 & mean(ev_rwd) < 0.95, label = paste(sc, "RWD events"))
  }
})

test_that("Sc4 with delta_rwd = 0 equals Sc1; Sc4 shifts only outside R; Sc5 everywhere", {
  for (r in 1:5) {
    d1 <- gen_data_id("Sc1", "gaussian", seed = r, n_mc = 0)
    d4 <- gen_data_id("Sc4_X5_d0", "gaussian", seed = r, n_mc = 0)
    expect_equal(d1$rct_trt, d4$rct_trt)
    expect_equal(d1$rct_ctrl, d4$rct_ctrl)
    expect_equal(d1$rwd_ctrl, d4$rwd_ctrl)
    expect_equal(d1$truth, d4$truth)
    d4s <- gen_data_id("Sc4_X5_d2", "gaussian", seed = r, n_mc = 0)
    inR <- d4s$rwd_ctrl$X5 > 2
    expect_equal(d4s$rwd_ctrl$y[inR], d1$rwd_ctrl$y[inR])
    expect_equal(d4s$rwd_ctrl$y[!inR], d1$rwd_ctrl$y[!inR] + 2)
    expect_equal(d4s$rct_ctrl, d1$rct_ctrl)   # RCT unchanged
    expect_equal(d4s$mu0_rct, d1$mu0_rct)
    expect_equal(d4s$mu0_rwd_at_rct - d4s$mu0_rct, 2 * (!d4s$region_rct))
    d5 <- gen_data_id("Sc5_d2", "gaussian", seed = r, n_mc = 0)
    expect_equal(d5$rwd_ctrl$y, d1$rwd_ctrl$y + 2)
    expect_true(all(!d5$region_rct))
    # secondary region
    d47 <- gen_data_id("Sc4_X7_d1", "gaussian", seed = r, n_mc = 0)
    inR7 <- d47$rwd_ctrl$X7 <= 2
    expect_equal(d47$rwd_ctrl$y[!inR7], d1$rwd_ctrl$y[!inR7] + 1)
    # survival Sc4 equals Sc1 at delta 0 as well
    s1 <- gen_data_id("Sc1", "survival", seed = r, n_mc = 0)
    s4 <- gen_data_id("Sc4_X5_d0", "survival", seed = r, n_mc = 0)
    expect_equal(s1$rwd_ctrl, s4$rwd_ctrl)
  }
  # region shares: R = {X5 > 2} holds about half of the RCT population
  sh <- vapply(seeds, function(r) mean(gen_data_id("Sc4_X5_d1", "gaussian", seed = r, n_mc = 0)$region_rct), 0)
  expect_true(abs(mean(sh) - 0.5) < 0.03)
})

test_that("reproducibility: same seed gives identical data", {
  a <- gen_data_id("Sc3_rho0.5", "survival", seed = 99, n_mc = 1e3)
  b <- gen_data_id("Sc3_rho0.5", "survival", seed = 99, n_mc = 1e3)
  expect_equal(a, b)
})

cat("test_dgp.R: all checks passed\n")
