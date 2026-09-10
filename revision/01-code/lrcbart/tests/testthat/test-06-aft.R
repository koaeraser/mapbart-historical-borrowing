# Test 6 (spec test 7): AFT with no censoring equals the Gaussian fit.
test_that("AFT without censoring reproduces the Gaussian fit exactly", {
  set.seed(9)
  n <- 150; x <- matrix(rnorm(n * 3), n, 3); src <- rep(1:2, c(50, 100))
  y <- 0.5 * x[, 1] + rnorm(n, 0, 0.5)
  fg <- lrc_bart(y, x, src, s0_sq = 0.01, n_burn = 200, n_draw = 200, seed = 10)
  fa <- lrc_bart_aft(exp(y), rep(1, n), x, src, s0_sq = 0.01, n_burn = 200, n_draw = 200, seed = 10)
  expect_equal(fa$sigma1_sq, fg$sigma1_sq)
  expect_equal(fa$f_train, fg$f_train)
  expect_equal(fa$g_train, fg$g_train)
  expect_equal(fa$tau0_sq, fg$tau0_sq)
  cat("\n[T6] AFT (no censoring) and Gaussian draws identical: TRUE\n")
})

test_that("AFT with censoring imputes above the censoring bound and recovers the scale", {
  set.seed(12)
  n <- 300; x <- matrix(rnorm(n * 3), n, 3); src <- rep(1:2, c(100, 200))
  logt <- 0.5 * x[, 1] + rnorm(n, 0, 0.5)
  cens <- log(rexp(n, 0.3))
  status <- as.integer(logt <= cens); time <- exp(pmin(logt, cens))
  fa <- lrc_bart_aft(time, status, x, src, s0_sq = 0.01, n_burn = 300, n_draw = 300, seed = 13)
  cat(sprintf("[T6b] censoring %.0f%%; sigma1 %.3f sigma2 %.3f (truth 0.5)\n",
              100 * mean(status == 0), sqrt(mean(fa$sigma1_sq)), sqrt(mean(fa$sigma2_sq))))
  expect_lt(abs(sqrt(mean(fa$sigma2_sq)) - 0.5), 0.1)
  mu <- colMeans(fa$f_train + (src == 1) * fa$g_train)
  expect_lt(sqrt(mean((mu - 0.5 * x[, 1])^2)), 0.25)
})
