# Test 2: two-leaf toy, one compatible leaf (residual 0) and one shifted
# (residual -3 sigma1), n1 = 50 per leaf.
test_that("two-leaf toy: P(spike) > 0.9 compatible, < 0.1 shifted (analytic block)", {
  sig1 <- 1.5; sig1sq <- sig1^2; tau0sq <- 0.01; tau1sq <- 1.25; w <- 0.8; n1 <- 50
  set.seed(2)
  pA <- mean(g_leaf_draw_cpp(n1, 0 * n1, sig1sq, tau0sq, tau1sq, w, 20000)$z)
  pB <- mean(g_leaf_draw_cpp(n1, -3 * sig1 * n1, sig1sq, tau0sq, tau1sq, w, 20000)$z)
  cat(sprintf("\n[T2 analytic] P(spike | A) = %.3f, P(spike | B) = %.4f\n", pA, pB))
  expect_gt(pA, 0.9)
  expect_lt(pB, 0.1)
})

test_that("two-leaf toy: full sampler with H_g = 1 separates the leaves", {
  set.seed(3)
  sig1 <- 1.5; sig2 <- 1.5
  n1 <- 100; n2 <- 200
  x1 <- c(runif(n1 / 2, -1, 0), runif(n1 / 2, 0, 1)); x2 <- runif(n2, -1, 1)
  y1 <- 1 + ifelse(x1 > 0, -3 * sig1, 0) + rnorm(n1, 0, sig1)
  y2 <- 1 + rnorm(n2, 0, sig2)
  fit <- lrc_bart(c(y1, y2), matrix(c(x1, x2), ncol = 1), c(rep(1, n1), rep(2, n2)),
                  H_f = 50, H_g = 1, w = 0.8, s0_sq = 0.01, n_burn = 500, n_draw = 1000, seed = 4)
  # region averages over a grid (pointwise f wiggles at noise sd 1.5 are not the target)
  xA <- matrix(seq(-0.9, -0.1, by = 0.05), ncol = 1); xB <- matrix(seq(0.1, 0.9, by = 0.05), ncol = 1)
  bA <- mean(lrc_borrowing_map(fit, xA, type = "all_spike")$mean); bB <- mean(lrc_borrowing_map(fit, xB, type = "all_spike")$mean)
  gA <- mean(lrc_discrepancy(fit, xA)); gB <- mean(lrc_discrepancy(fit, xB))
  cat(sprintf("[T2 sampler] mean B in A = %.3f, in B = %.4f; mean g in A = %.3f, in B = %.3f (truth 0, -4.5; RCT sample means %.2f, %.2f)\n",
              bA, bB, gA, gB, mean(y1[x1 < 0]) - mean(y2), mean(y1[x1 > 0]) - mean(y2)))
  expect_gt(bA, 0.8)
  expect_lt(bB, 0.1)
  expect_lt(abs(gB + 4.5), 0.6)
  expect_lt(abs(gA), 0.3)
})
