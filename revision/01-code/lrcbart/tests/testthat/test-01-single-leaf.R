# Test 1: single g leaf, f fixed. Sampler P(z = 1) and E[theta] match T1.
test_that("single g-leaf block draw matches the T1 formulas over a grid", {
  sig1sq <- 1.5^2; tau0sq <- 0.01; tau1sq <- 1.25; w <- 0.8
  t1 <- function(ybar, n1) {
    v1 <- sig1sq / n1
    lo <- log((1 - w) / w) + 0.5 * log((tau0sq + v1) / (tau1sq + v1)) +
      0.5 * (1 / (tau0sq + v1) - 1 / (tau1sq + v1)) * ybar^2
    p <- plogis(lo)                       # P(z = 0 | ybar), the slab probability
    th1 <- tau1sq * ybar / (tau1sq + v1)  # slab estimate
    th0 <- tau0sq * ybar / (tau0sq + v1)  # spike estimate
    c(pspike = 1 - p, mean = p * th1 + (1 - p) * th0)
  }
  set.seed(11)
  grid <- expand.grid(ybar = c(-2, -1, -0.5, -0.2, 0, 0.3, 1, 3), n1 = c(5, 20, 50, 200))
  ndraw <- 40000
  worst_p <- 0; worst_m <- 0
  for (r in seq_len(nrow(grid))) {
    yb <- grid$ybar[r]; n1 <- grid$n1[r]
    ref <- t1(yb, n1)
    p_cpp <- g_leaf_pspike_cpp(n1, yb * n1, sig1sq, tau0sq, tau1sq, w)
    expect_equal(p_cpp, unname(ref["pspike"]), tolerance = 1e-10)
    d <- g_leaf_draw_cpp(n1, yb * n1, sig1sq, tau0sq, tau1sq, w, ndraw)
    se_p <- sqrt(ref["pspike"] * (1 - ref["pspike"]) / ndraw) + 1e-6
    se_m <- sd(d$theta) / sqrt(ndraw) + 1e-6
    worst_p <- max(worst_p, abs(mean(d$z) - ref["pspike"]) / se_p)
    worst_m <- max(worst_m, abs(mean(d$theta) - ref["mean"]) / se_m)
    expect_lt(abs(mean(d$z) - ref["pspike"]), 4.5 * se_p)
    expect_lt(abs(mean(d$theta) - ref["mean"]), 4.5 * se_m)
  }
  cat(sprintf("\n[T1] worst |P(z) error| = %.2f SE, worst |E(theta) error| = %.2f SE over %d grid points\n",
              worst_p, worst_m, nrow(grid)))
})
