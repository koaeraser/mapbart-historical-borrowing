# Test 5 (spec test 6): tau0^2, tau1^2 and w updates recover known values.
test_that("update_tau0 / update_tau1 / update_w recover known scales", {
  set.seed(8)
  L <- 2000; K0 <- 1500
  z <- c(rep(1L, K0), rep(0L, L - K0))
  theta <- c(rnorm(K0, 0, sqrt(0.05)), rnorm(L - K0, 0, sqrt(2)))
  t0 <- update_tau0_cpp(theta, z, 3, 0.01, 10, 4000)
  t1 <- update_tau1_cpp(theta, z, 3, 1, 0.05, 4000)
  ww <- update_w_cpp(z, 4, 1, 4000)
  cat(sprintf("\n[T5] tau0^2: truth 0.050, posterior mean %.4f (S0/K0 = %.4f); tau1^2: truth 2.00, posterior mean %.3f (S1/K1 = %.3f); w: truth 0.75, posterior mean %.3f\n",
              mean(t0), sum(theta[z == 1]^2) / K0, mean(t1), sum(theta[z == 0]^2) / (L - K0), mean(ww)))
  expect_lt(abs(mean(t0) / 0.05 - 1), 0.1)
  expect_lt(abs(mean(t1) / 2 - 1), 0.1)
  expect_lt(abs(mean(ww) - (4 + K0) / (5 + L)), 0.01)
  # truncation is respected
  expect_true(all(t0 < 10)); expect_true(all(t1 > 0.05))
  t0b <- update_tau0_cpp(theta, z, 3, 0.01, 0.02, 4000)   # truncate well below the data scale
  expect_true(all(t0b < 0.02))
  # exact conjugate check with no data: prior mean nu0 s0 / (nu0 - 2)
  pr <- update_tau0_cpp(numeric(0), integer(0), 5, 0.1, 1e6, 200000)
  expect_lt(abs(mean(pr) / (5 * 0.1 / 3) - 1), 0.05)
  # truncated normal sampler: mean of N(0,1) truncated at a = 2 is phi(2)/(1-Phi(2))
  r <- rtnorm_lower_cpp(0, 1, 2, 100000)
  expect_true(all(r > 2))
  expect_lt(abs(mean(r) - dnorm(2) / pnorm(2, lower.tail = FALSE)), 0.01)
})
