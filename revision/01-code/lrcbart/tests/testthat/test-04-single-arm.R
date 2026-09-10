# Test 4 (spec test 5): single-arm w = 1 reproduces standard BART on the RWD.
test_that("single-arm w = 1: f equals plain BART on the RWD and matches BART::wbart", {
  skip_if_not_installed("BART")
  set.seed(6)
  n <- 200; p <- 5
  x <- matrix(rnorm(n * p), n, p)
  ftrue <- function(x) 2 * sin(x[, 1]) + x[, 2]^2 - x[, 3]
  y <- ftrue(x) + rnorm(n)
  xt <- matrix(rnorm(100 * p), 100, p)
  # single-arm LRC-BART (no RCT controls), w = 1
  fit <- lrc_bart(y, x, rep(2, n), w = 1, s0_sq = 0.01, n_burn = 1000, n_draw = 1000, seed = 7)
  expect_true(fit$single_arm)
  # in single-arm mode g is a prior draw: mean 0, pointwise variance H_g tau0^2
  g <- lrc_discrepancy(fit, xt)
  expect_lt(abs(mean(g)), 0.05)
  expect_equal(mean(apply(g, 2, var)), mean(fit$H_g * fit$tau0_sq), tolerance = 0.3)
  # f draws equal a plain bart_fit on the same data and seed (same RNG stream
  # up to the g prior draws, so compare means within Monte Carlo error)
  pb <- bart_fit(y, x, n_burn = 1000, n_draw = 1000, seed = 7)
  m_lrc <- colMeans(predict(fit, xt, arm = "rwd_control"))
  m_pb <- colMeans(predict(pb, xt))
  # wbart, 50 trees, same defaults (k = 2, nu = 3, q = 0.9, 100 cutpoints)
  wb <- BART::wbart(x, y, xt, ntree = 50, ndpost = 1000, nskip = 1000, printevery = 1e6)
  m_wb <- wb$yhat.test.mean
  sd_post <- mean(apply(predict(pb, xt), 2, sd))
  cat(sprintf("\n[T4] rmse vs truth: lrc single-arm %.3f, bart_fit %.3f, wbart %.3f; cor(lrc, wbart) = %.4f; mean|diff| = %.3f (posterior sd %.3f); sigma: lrc %.3f wbart %.3f\n",
              sqrt(mean((m_lrc - ftrue(xt))^2)), sqrt(mean((m_pb - ftrue(xt))^2)),
              sqrt(mean((m_wb - ftrue(xt))^2)), cor(m_lrc, m_wb), mean(abs(m_lrc - m_wb)), sd_post,
              sqrt(mean(fit$sigma2_sq)), mean(wb$sigma)))
  expect_gt(cor(m_lrc, m_wb), 0.95)
  expect_lt(mean(abs(m_lrc - m_wb)), 0.5 * sd_post)
  expect_lt(mean(abs(m_lrc - m_pb)), 0.5 * sd_post)
  expect_lt(abs(sqrt(mean(fit$sigma2_sq)) - mean(wb$sigma)), 0.1)
  # single-arm sigma1 is tied to sigma2
  expect_equal(fit$sigma1_sq, fit$sigma2_sq)
})
