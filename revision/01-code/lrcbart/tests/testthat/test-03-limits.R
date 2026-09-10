# Test 3 (spec test 4): w = 0 and w = 1 limits of the g leaf factor.
test_that("w = 0 gives the free-offset leaf factor; w = 1 gives C-BART; tau0 -> 0 gives BART-CP", {
  sig1sq <- 2; n1 <- 30; s1 <- 12; tau0sq <- 0.02; tau1sq <- 3
  free <- f_leaf_logmarg_cpp(n1, 0, s1, 0, sig1sq, 1, tau1sq)
  expect_equal(g_leaf_logmarg_cpp(n1, s1, sig1sq, tau0sq, tau1sq, 0), free, tolerance = 1e-10)
  cb <- f_leaf_logmarg_cpp(n1, 0, s1, 0, sig1sq, 1, tau0sq)
  expect_equal(g_leaf_logmarg_cpp(n1, s1, sig1sq, tau0sq, tau1sq, 1), cb, tolerance = 1e-10)
  expect_equal(g_leaf_logmarg_cpp(n1, s1, sig1sq, 1e-14, tau1sq, 1), 0, tolerance = 1e-8)
  # tau1 -> infinity under w = 0: the leaf factor tends to that of a flat prior,
  # -0.5 log(1 + tau1^2 n1 / sig1^2) + ..., i.e. log M -> -Inf like -0.5 log tau1^2
  # while the difference between two ybar values converges to the likelihood ratio
  big <- 1e8
  d <- g_leaf_logmarg_cpp(n1, s1, sig1sq, tau0sq, big, 0) - g_leaf_logmarg_cpp(n1, 0, sig1sq, tau0sq, big, 0)
  expect_equal(d, (s1 / n1)^2 / (2 * sig1sq / n1), tolerance = 1e-5)
  # P(z) at the limits
  expect_equal(g_leaf_pspike_cpp(n1, s1, sig1sq, tau0sq, tau1sq, 1), 1)
  expect_equal(g_leaf_pspike_cpp(n1, s1, sig1sq, tau0sq, tau1sq, 0), 0)
})

test_that("w = 1 and w = 0 are reachable in the sampler", {
  set.seed(5)
  n <- 120; x <- matrix(rnorm(n * 2), n, 2); src <- rep(1:2, each = n / 2)
  y <- x[, 1] + rnorm(n)
  f1 <- lrc_bart(y, x, src, w = 1, s0_sq = 0.01, n_burn = 50, n_draw = 50, seed = 1)
  expect_true(all(f1$K0 == f1$L_g)); expect_true(all(f1$w == 1))
  f0 <- lrc_bart(y, x, src, w = 0, s0_sq = 0.01, n_burn = 50, n_draw = 50, seed = 1)
  expect_true(all(f0$K0 == 0)); expect_true(all(f0$w == 0))
  # w = 0 borrowing map is identically 0, w = 1 identically 1
  expect_equal(unname(lrc_borrowing_map(f0, x[1:3, ], type = "all_spike")$mean), c(0, 0, 0))
  expect_equal(unname(lrc_borrowing_map(f1, x[1:3, ], type = "all_spike")$mean), c(1, 1, 1))
})
