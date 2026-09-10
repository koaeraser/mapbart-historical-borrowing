# Test 8 (Phase 2b options): fractional ESS target, slab multiplier, g leaf
# minimum, extra g sweeps and a warm start all run and do what they say.
test_that("fractional N_target, tau1_mult, n_min_g, g_sweeps and g_init behave", {
  set.seed(8)
  n1 <- 60; n2 <- 150; p <- 4
  x1 <- matrix(rnorm(n1 * p), n1, p); x2 <- matrix(rnorm(n2 * p), n2, p)
  f <- function(X) X[, 1] + 0.5 * X[, 2]^2
  y1 <- f(x1) - 1.5 * (x1[, 1] < 0) + rnorm(n1); y2 <- f(x2) + rnorm(n2, 0, 1.2)
  cal <- lrc_ess_calibrate(y2, x2, rbind(x1, x2), N_target = c("frac", 0.9), y_rct_ctrl = y1, x_rct_ctrl = x1,
                           n_burn = 200, n_draw = 400, seed = 3, warn = FALSE)
  expect_false(cal$ess_capped)
  expect_equal(cal$N_target, 0.9 * cal$ceiling)
  expect_identical(as.character(cal$N_target_spec[[1]]), "frac")
  cal2 <- lrc_ess_calibrate(y2, x2, rbind(x1, x2), N_target = list("frac", 0.5), y_rct_ctrl = y1, x_rct_ctrl = x1,
                            n_burn = 200, n_draw = 400, seed = 3, warn = FALSE, fit_stage1 = cal$fit_stage1)
  expect_equal(cal2$N_target, 0.5 * cal2$ceiling)
  expect_gte(cal2$s0_sq, cal$s0_sq)  # a lower target needs a tighter spike or the same grid value
  y <- c(y1, y2); x <- rbind(x1, x2); src <- rep(1:2, c(n1, n2))
  f1 <- lrc_bart(y, x, src, s0_sq = cal$s0_sq, calib = cal, n_burn = 200, n_draw = 200, seed = 1, tau1_mult = 4)
  f0 <- lrc_bart(y, x, src, s0_sq = cal$s0_sq, calib = cal, n_burn = 200, n_draw = 200, seed = 1)
  expect_equal(f1$hyper$tau1_sq, 4 * f0$hyper$tau1_sq)
  expect_equal(f0$hyper$H_g, 5L); expect_equal(f0$hyper$a_w, 1); expect_equal(f0$hyper$b_w, 1)
  f2 <- lrc_bart(y, x, src, s0_sq = cal$s0_sq, calib = cal, n_burn = 200, n_draw = 200, seed = 1, n_min_g = 1, g_sweeps = 3)
  expect_equal(f2$hyper$n_min_g, 1L)
  expect_equal(length(f2$draws$g), 200L)  # leaf rows carry NA cut values by design
  # warm start: a greedy tree on the true offset, loaded into tree 1 of g
  off <- -1.5 * (x1[, 1] < 0)
  gi <- lrc_g_init_greedy(x1, off, H_g = 5, maxdepth = 1)
  expect_equal(length(gi$roots), 5)
  expect_equal(unname(gi$nodes[1, "var"]), 1)
  f3 <- lrc_bart(y, x, src, s0_sq = cal$s0_sq, calib = cal, n_burn = 0, n_draw = 1, seed = 1, g_init = gi)
  g_first <- lrc_discrepancy(f3, x1)
  # after one iteration from the warm start the sign pattern of the offset is present
  expect_lt(mean(g_first[1, x1[, 1] < 0]), mean(g_first[1, x1[, 1] >= 0]))
})
