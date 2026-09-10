# test_cahb.R: checks for the CAHB comparator (R/cahb.R).
#   (i)   Sc1 (compatible RWD): CAHB borrows (map, realized ESS), its ATE SD
#         is below its own no-borrowing SD, and (as asked) below BART-NP's.
#         The last check is expected to fail: CAHB's kernel smoother leaves a
#         residual variance of 4 to 5 against BART's 2.3, so its posterior
#         SD is larger than BART-NP's even with borrowing. The result is
#         printed either way.
#   (ii)  Sc5_d2 (global shift of 2): CAHB borrows much less than in Sc1.
#   (iii) Draw dimensions, no NA, map in [0, 1], lengths.
# Then per-fit timing with time_methods().
# Usage: Rscript tests/test_cahb.R
suppressPackageStartupMessages(library(testthat))
CODE_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), ".."))
source(file.path(CODE_ROOT, "R", "harness.R")); source_all(CODE_ROOT)
source(file.path(CODE_ROOT, "R", "cahb.R"))

seeds <- c(3001, 3002)   # outside the study's data-seed range 2027..2526
fits <- list()
for (sc in c("Sc1", "Sc5_d2")) for (s in seeds) {
  dat <- gen_data_id(sc, "gaussian", seed = s, n_mc = 0)
  set.seed(s); f_c <- fit_cahb(dat, "gaussian")
  set.seed(s); f_nb <- fit_cahb(dat, "gaussian", lam2_coef = 0)   # lambda_2 = 0: no borrowing
  set.seed(s); f_b <- fit_bart_np(dat, "gaussian")
  fits[[paste(sc, s)]] <- list(dat = dat, cahb = f_c, nb = f_nb, bart = f_b,
                               m_c = compute_metrics(f_c, dat), m_nb = compute_metrics(f_nb, dat),
                               m_b = compute_metrics(f_b, dat))
}
show <- function(k) {
  x <- fits[[k]]
  cat(sprintf("%-12s CAHB: ate %.3f sd %.3f map_all %.3f ess %.1f surf %.3f | no-borrow sd %.3f | BART-NP: ate %.3f sd %.3f | kernel vars %s mult %.2f phi2_hat %.2f lam1 %.2f lam2 %.1f\n",
              k, x$m_c$ate_pm, x$m_c$ate_sd, x$m_c$map_all, x$m_c$ess_realized, x$m_c$surf_rmse,
              x$m_nb$ate_sd, x$m_b$ate_pm, x$m_b$ate_sd,
              paste0("X", x$cahb$extra$vars, collapse = ","), x$cahb$extra$mult_c,
              x$cahb$extra$phi2_hat, x$cahb$extra$lam1, x$cahb$extra$lam2))
}
for (k in names(fits)) show(k)

with_reporter("summary", {
test_that("(i) Sc1: CAHB borrows and borrowing reduces its own ATE SD", {
  for (s in seeds) {
    x <- fits[[paste("Sc1", s)]]
    expect_gt(x$m_c$map_all, 0.2)
    expect_gt(x$m_c$ess_realized, 30)
    expect_lt(x$m_c$ate_sd, x$m_nb$ate_sd)
  }
})

test_that("(i-c) Sc1: CAHB ATE SD below BART-NP's (expected to fail, see header)", {
  for (s in seeds) {
    x <- fits[[paste("Sc1", s)]]
    expect_lt(x$m_c$ate_sd, x$m_b$ate_sd)
  }
})

test_that("(ii) Sc5_d2: CAHB borrows much less than under Sc1", {
  for (s in seeds) {
    x1 <- fits[[paste("Sc1", s)]]; x5 <- fits[[paste("Sc5_d2", s)]]
    expect_lt(x5$m_c$map_all, 0.1)
    expect_lt(x5$m_c$ess_realized, 0.2 * x1$m_c$ess_realized)
  }
})

test_that("(iii) dimensions, no NA, map range", {
  for (k in names(fits)) {
    x <- fits[[k]]; f <- x$cahb; dat <- x$dat
    expect_equal(dim(f$draws), c(N_DRAW, 3L))
    expect_named(f$draws, c("ate", "trt", "ctrl"))
    expect_false(anyNA(f$draws))
    expect_length(f$surface, nrow(dat$x_rct))
    expect_false(anyNA(f$surface))
    expect_length(f$map, sum(dat$z_rct == 0))
    expect_true(all(f$map >= 0 & f$map < 1))
    expect_true(is.finite(f$ess$realized))
    prov <- setdiff(grep("^(ate|trt|ctrl)_|^surf_rmse$|^map_all$|^ess_realized$", names(x$m_c), value = TRUE),
                    c("trt_reject", "trt_reject90", "ctrl_reject", "ctrl_reject90"))   # NA by construction (no null for arm means)
    expect_true(all(is.finite(unlist(x$m_c[prov]))))   # columns CAHB provides (others are NA by design)
  }
})
})   # with_reporter: failures are reported, not fatal

cat("\nTiming (mean seconds per fit, 3 replicates of Sc1, one core):\n")
tm <- time_methods("Sc1", "gaussian", c("CAHB", "CAHB-p10", "CAHB-lit", "BART-NP"), n_rep = 3, seed = 4000)
print(tm[, c("method", "mean_seconds", "sd_seconds")], row.names = FALSE)
