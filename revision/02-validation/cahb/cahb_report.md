# CAHB comparator: implementation, tests and paired results

Prepared 2026-09-09 for the LRC-BART revision. Method: covariate-adaptive historical control borrowing (CAHB) of Jin, Kim, Scheffler and Jiang (2023), Statistics in Medicine 42(29):5338-5352, doi 10.1002/sim.9913. Code: `01-code/R/cahb.R`; tests: `01-code/tests/test_cahb.R`; run and comparison: this directory (`run_cahb.R`, `compare_cahb.R`, `replay_extras.R`, `results/`, `summary_table.csv`, `paired.csv`, `extras_cahb.csv`, `summary_table_cahb_only.csv`, `wall_time.txt`). No existing file in `01-code/R` was modified; `cahb.R` registers its methods when sourced after `comparators.R`, and the run script sources it explicitly.

## 1. Source and what CAHB is

We obtained the published paper and its supplement from the first author's site (jinhuaqing.github.io/files/2023_SIM_CAHB.pdf and 2023_SIM_CAHB_supp.pdf; the PMC copy PMC10919261 and the Wiley page are behind a CAPTCHA and a 403 respectively), and the authors' code from github.com/JINhuaqing/HistTrial (simulation scripts, not a package; nothing on CRAN). We read the method sections, the supplement's closed-form posterior, and the code's `utils.R`, `simplex.R` and `realData/linearSimuRealData*.R` in full.

CAHB is a sequential trial design. Its analysis component places a local prior on the current control mean function, mu_0(X) ~ N(theta_0(X), 1/tau(X)), where theta_0 is the historical control mean function supplied as a known input and tau(X) is a covariate-dependent precision with a half-normal prior of scale gamma on tau^{1/2}; phi_0^2 and phi_1^2 have Inv-Gamma(0.01, 0.01) priors and mu_1 a flat prior. Everything is estimated by kernel-local posterior modes iterated to convergence (Algorithm 1): mu_0(X_j) is a kernel-weighted combination of the control outcomes and the pseudo-outcomes theta_0(X_i) weighted by tau_i; tau~(X_j) is the local ratio of the kernel mass to the kernel-weighted squared discrepancy {mu_0(X_i) - theta_0(X_i)}^2 plus 1/gamma^2; the vector of tau~ over the trial subjects is thresholded at lambda_1 (Algorithm 2: the 10% quantile of tau~ when theta_0 is replaced by the no-borrowing fit of the current controls) and projected onto the L1 ball of radius lambda_2 = 300 log n; phi_0 is the posterior mode of the residual variance. The plug-in posterior of mu_0(X) is normal with variance 1/A(X), A(X) = sum_i (1-Z_i) w_ij {1/phi_0^2 + tau_i}, and the information ratio R(X) = A(X) phi_0^2 / sum_i (1-Z_i) w_ij (supplement A.2) defines the covariate-dependent effective control sample size, R(X) times the local control count. The design part (a kernel biased-coin allocation that assigns more subjects to the experimental arm where R(X) is large) has no role in a completed RCT with fixed allocation and is not implemented. The published constants are gamma = sqrt(3), lambda_2 = 300 log(n), the 10% quantile for lambda_1, a Gaussian product kernel with Scott's rule-of-thumb bandwidths for continuous covariates (0.1 for binary ones), and convergence at 1e-5 within 100 iterations. The paper's own settings have four covariates (two binary), a residual variance of 0.11 (24-month hip BMD), N = 200 current subjects and a historical model that is linear within the four binary-defined subgroups.

## 2. What we implemented, and every choice the paper leaves open

The primary method "CAHB" is the analysis model above with the following choices, all recorded in the header of `cahb.R`. Two variants run alongside: "CAHB-p10" (all ten covariates in the kernel, otherwise identical) and "CAHB-lit" (every published constant taken literally). The choices were settled on three diagnostic replicates of Sc1, Sc5_d2 and Sc4_X5_d2 (seeds 2027, 2050, 2101; the study seeds are 2026 + r and the test seeds 3001, 3002) by whether the estimator operates at all: local sample size, the phi_0 estimate, and whether the tau switch responds to a global shift. They were not tuned to the study metrics.

1. Historical mean function theta_0. The paper fits a subgroup-wise linear model to the historical patients. We use the same kernel regression the method uses for mu_0 (Nadaraya-Watson on the kernel covariates, same bandwidth), fitted to the 300 RWD controls and evaluated at the RCT profiles. A sharper historical estimate (BART on the RWD) was tried and rejected: paired with the kernel-smoothed mu_0, the discrepancy mu_0 - theta_0 then reflects the smoother's bias rather than the shift, tau is small everywhere and nothing is borrowed even under Sc1. As in the paper, theta_0 enters as a known function and its variance is ignored (the paper states this as a limitation in its Section 6).

2. Kernel covariates. The paper matches on four literature-chosen prognostic covariates. We use the four covariates with the largest variable-inclusion counts in a BART fit to the RWD controls (the harness's BART_OPTS; BART is used only for this ranking). Over the 1000 study fits the kernel always contained X1, X5 and X6, plus X3 (32%), X4 (32%), X2 (14%) or one of X7 to X10 (22%); `extras_cahb.csv` lists the set per replicate.

3. Bandwidth. Scott's rule h_k = s_k n^{-1/(p+4)} with p = 10 and the n = 100 current controls gives h = 0.72 s_k, at which the product kernel localizes to single subjects: the median local effective control sample size sum_i w_ij is about 1.0, of which 1 is the subject itself. The estimator then interpolates, phi_0 collapses to about 0.1, profiles without a neighbour get a near-infinite posterior variance, and the tau threshold cannot separate a shift from noise, since separation needs sum_i w_ij d^2 much larger than 1/gamma^2. Leave-one-out cross-validation of the no-borrowing smoother does not repair this: in four or ten dimensions it selects 1 to 1.5 times Scott's rule (local sample 1 to 5) because Nadaraya-Watson is poor at every bandwidth for this surface, and at that bandwidth CAHB borrowed from the shifted RWD in Sc5_d2 on two of three diagnostic replicates (ATE estimates of -0.04 and 0.12 for a truth of 1). We therefore multiply Scott's rule by the constant at which the median local effective control sample size over the RCT control profiles equals 20, found by bisection (2.5 in four dimensions, 3.4 in ten; the smallest local sample in a replicate is then about 1.2). The kernel is normalized to K(0) = 1 as in the paper's Section 3, so sum_i w_ij is a local sample size; the authors' code uses an unnormalized Gaussian density, which rescales the posterior variance by (2 pi)^{-p/2} |H|^{-1/2} (about 0.7 in their real-data setting) and alters tau~ through the 1/gamma^2 term. One H serves mu_0, tau and mu_1, as in the paper.

4. gamma and lambda_2 are on the outcome scale. The published gamma = sqrt(3) and lambda_2 = 300 log n go with the paper's residual variance 0.11; the authors' own linear-model script (`linearSimu.R`, phi_0 = 3) uses gamma^2 = 0.01 and lambda_2 = 30 log n. The dimensionless quantities are gamma^2 phi_0^2, the prior cap on the borrowed precision relative to the data precision, and lambda_2 phi_0^2 / n, the projection cap on the mean borrowed-to-data precision ratio. We keep the published dimensionless values, gamma^2 = 3 x 0.11 / phi_hat^2 and lambda_2 = 300 log(n) x 0.11 / phi_hat^2, with phi_hat^2 the no-borrowing residual variance of the RCT controls under the method's own smoother (the reference model of Algorithm 2). In the study phi_hat^2 averaged 5.1 (the true residual variance is 2.25; the difference is the structure the smoother leaves in the residual), giving gamma^2 = 0.066, lambda_2 = 37.5 and lambda_1 = 0.45 on average (`extras_cahb.csv`).

5. Posterior draws of the averaged estimand. The authors' code draws mu_0(X_j) and mu_1(X_j) independently across profiles from their marginal posteriors and averages; the resulting posterior of the profile-averaged difference is not calibrated, and the paper calibrates its test cut-off c by simulation instead. We keep every marginal posterior of the paper (mean mu_0(X_j), variance 1/A(X_j)) and induce the dependence across profiles through the smoother's own weights: with a_ij = (1-Z_i) w_ij / phi_0^2 / A(X_j) and standard normal epsilon_i, one draw is mu_0(X_j) + A(X_j)^{-1/2} sum_i a_ij epsilon_i / (sum_i a_ij^2)^{1/2}, and likewise for mu_1. Point estimates, the surface and the map do not depend on this; the posterior SD, coverage and power do.

6. Treated arm and standardization: mu_1 from the RCT treated subjects only (flat prior), estimands standardized over all n RCT profiles as for every other method; the paper's point estimate N^{-1} sum_i {mu_1(X_i) - mu_0(X_i)} is the same quantity.

7. tau is evaluated at all n current-trial profiles and the L1 projection acts on that length-n vector, as in the authors' `mOptTau`.

8. Harness outputs: map = 1 - 1/R(X) at the RCT control profiles (the historical share of the local posterior precision, in [0, 1)); ess = list(realized = sum over RCT controls of R(X_i) - 1, the CECSS excess over the actual control count; no prior ESS, since CAHB has none); extra carries tau, R, the constants, the kernel covariates and counts, the local sample sizes and the iteration count.

CAHB-lit takes the published values literally: theta_0 by linear regression on the RWD (the paper's historical model), all ten covariates, Scott's rule as is, gamma^2 = 3, lambda_2 = 300 log n, independent per-profile draws. Gaussian outcome only; the paper treats a continuous outcome and we did not attempt the survival study.

## 3. Tests and timing

`Rscript tests/test_cahb.R` (seeds 3001, 3002). Test (i), Sc1: CAHB borrows (map_all 0.53 and 0.53, realized ESS 116 and 128) and its ATE SD with borrowing (0.42, 0.40) is below its own no-borrowing SD (0.50, 0.47): pass. The check asked for in the task, ATE SD below BART-NP's, fails on both replicates (0.42 vs 0.19, 0.40 vs 0.18); the test is kept in the file and marked as expected to fail. The reason is that CAHB's kernel smoother leaves a residual variance of 4 to 5 against BART's 2.3, and its posterior SD scales with that residual variance. Test (ii), Sc5_d2: map_all 0.000 and realized ESS 0 on both replicates against 116 and 128 under Sc1: pass. Test (iii): draws 1000 x 3 named ate/trt/ctrl, no NA in draws, surface (length 300) or map (length 104, all in [0, 1)), all metric columns that CAHB provides finite: pass. Timing by `time_methods()` on three Sc1 replicates, one core: CAHB 0.61 to 0.67 s per fit (0.5 s of it is the BART ranking), CAHB-p10 0.05 to 0.08 s, CAHB-lit 0.02 s, BART-NP 0.73 to 0.84 s. In the study the mean was 0.63 to 0.81 s for CAHB (`summary_table.csv`, column seconds).

## 4. Run

`OMP_NUM_THREADS=1 Rscript run_cahb.R`: Gaussian, scenarios Sc1, Sc4_X5_d1, Sc4_X5_d2, Sc5_d1, Sc5_d2, methods CAHB, CAHB-p10, CAHB-lit, R = 200, seed 2026, 12 workers; 3000 fits in 1.4 minutes wall, no fit errors. `compare_cahb.R` reads replicates 1 to 200 of `../full/gaussian/per_replicate.csv` for LRC-BART-100, BART-PP, BART-CP, BART-NP, PSCL and MAP-100, checks that ate_truth is identical per (scenario, replicate) in both files (1000 pairs, max difference below 1e-8), and writes `summary_table.csv` and `paired.csv`. `replay_extras.R` re-fits CAHB on the 1000 study replicates with the harness seeds and reproduces the stored ESS and ATE to 1e-12, so the diagnostics in `extras_cahb.csv` describe the stored run.

## 5. Results

The tables below are generated from `summary_table.csv` and `paired.csv`. Sc4 uses R = {X5 > 2}; in Sc5 the shift applies everywhere, so the R columns are empty.

### Table A. Paired summary, replicates 1 to 200 (Gaussian, seed 2026)

ATE bias, mean posterior SD, MC RMSE of the posterior mean, 95% coverage, power (95% CrI excludes 0), control-surface RMSE at the RCT profiles overall, in R and in R^c, mean borrowing map in R and R^c at the RCT control profiles, realized ESS (CAHB: CECSS excess; LRC-BART: posterior realized ESS), mean seconds per fit.

| scenario | method | bias | sd | rmse | cov | pow | surf_rmse | surf_rmse_R | surf_rmse_Rc | map_R | map_Rc | ess_realized | secs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | CAHB | -0.003 | 0.410 | 0.258 | 1.000 | 0.780 | 2.176 | 2.071 | 2.273 | 0.52 | 0.51 | 114.5 | 0.63 |
| Sc1 | CAHB-p10 | -0.018 | 0.542 | 0.262 | 1.000 | 0.335 | 2.408 | 2.298 | 2.509 | 0.50 | 0.50 | 103.1 | 0.06 |
| Sc1 | CAHB-lit | -0.087 | 12.226 | 0.950 | 1.000 | 0.065 | 4.627 | 4.826 | 3.160 | 0.01 | 0.01 | 0.6 | 0.02 |
| Sc1 | LRC-BART-100 | -0.009 | 0.171 | 0.175 | 0.935 | 1.000 | 0.757 | 0.740 | 0.771 | 0.92 | 0.92 | 52.3 | 2.08 |
| Sc1 | BART-PP | -0.014 | 0.207 | 0.206 | 0.940 | 1.000 | 0.773 | 0.754 | 0.789 |  |  |  | 1.49 |
| Sc1 | BART-CP | -0.006 | 0.145 | 0.148 | 0.960 | 1.000 | 0.751 | 0.731 | 0.768 |  |  |  | 1.55 |
| Sc1 | BART-NP | -0.016 | 0.202 | 0.223 | 0.895 | 1.000 | 1.160 | 1.146 | 1.169 |  |  |  | 1.05 |
| Sc1 | PSCL | -0.017 | 0.269 | 0.262 | 0.940 | 0.975 | 2.492 | 2.356 | 2.617 |  |  |  | 0.00 |
| Sc1 | MAP-100 | -0.022 | 0.290 | 0.266 | 0.950 | 0.945 | 2.492 | 2.356 | 2.617 |  |  |  | 0.90 |
| Sc4_X5_d1 | CAHB | -0.208 | 0.413 | 0.333 | 0.980 | 0.425 | 2.173 | 2.133 | 2.207 | 0.53 | 0.48 | 109.2 | 0.78 |
| Sc4_X5_d1 | CAHB-p10 | -0.260 | 0.545 | 0.380 | 1.000 | 0.125 | 2.423 | 2.384 | 2.456 | 0.50 | 0.49 | 101.2 | 0.13 |
| Sc4_X5_d1 | CAHB-lit | -0.088 | 12.207 | 0.948 | 1.000 | 0.065 | 4.621 | 4.820 | 3.158 | 0.01 | 0.01 | 0.6 | 0.02 |
| Sc4_X5_d1 | LRC-BART-100 | -0.064 | 0.194 | 0.227 | 0.870 | 1.000 | 0.820 | 0.778 | 0.856 | 0.64 | 0.53 | 21.6 | 2.41 |
| Sc4_X5_d1 | BART-PP | -0.023 | 0.209 | 0.203 | 0.960 | 1.000 | 0.803 | 0.776 | 0.826 |  |  |  | 1.69 |
| Sc4_X5_d1 | BART-CP | -0.352 | 0.146 | 0.382 | 0.315 | 0.995 | 0.897 | 0.750 | 1.021 |  |  |  | 1.70 |
| Sc4_X5_d1 | BART-NP | -0.016 | 0.202 | 0.223 | 0.895 | 1.000 | 1.160 | 1.146 | 1.169 |  |  |  | 1.15 |
| Sc4_X5_d1 | PSCL | -0.268 | 0.270 | 0.373 | 0.865 | 0.765 | 2.505 | 2.449 | 2.556 |  |  |  | 0.00 |
| Sc4_X5_d1 | MAP-100 | -0.254 | 0.303 | 0.381 | 0.890 | 0.715 | 2.505 | 2.442 | 2.561 |  |  |  | 1.02 |
| Sc4_X5_d2 | CAHB | -0.333 | 0.422 | 0.453 | 0.955 | 0.250 | 2.167 | 2.172 | 2.157 | 0.49 | 0.41 | 93.9 | 0.80 |
| Sc4_X5_d2 | CAHB-p10 | -0.265 | 0.585 | 0.533 | 0.975 | 0.170 | 2.421 | 2.386 | 2.448 | 0.31 | 0.31 | 62.5 | 0.13 |
| Sc4_X5_d2 | CAHB-lit | -0.089 | 12.187 | 0.947 | 1.000 | 0.065 | 4.616 | 4.815 | 3.155 | 0.00 | 0.00 | 0.5 | 0.02 |
| Sc4_X5_d2 | LRC-BART-100 | -0.020 | 0.188 | 0.206 | 0.905 | 1.000 | 0.844 | 0.820 | 0.862 | 0.66 | 0.05 | 0.9 | 2.21 |
| Sc4_X5_d2 | BART-PP | -0.029 | 0.210 | 0.203 | 0.950 | 1.000 | 0.809 | 0.780 | 0.832 |  |  |  | 1.57 |
| Sc4_X5_d2 | BART-CP | -0.697 | 0.149 | 0.714 | 0.005 | 0.525 | 1.239 | 0.772 | 1.573 |  |  |  | 1.57 |
| Sc4_X5_d2 | BART-NP | -0.016 | 0.202 | 0.223 | 0.895 | 1.000 | 1.160 | 1.146 | 1.169 |  |  |  | 1.06 |
| Sc4_X5_d2 | PSCL | -0.519 | 0.273 | 0.581 | 0.495 | 0.385 | 2.543 | 2.562 | 2.520 |  |  |  | 0.00 |
| Sc4_X5_d2 | MAP-100 | -0.410 | 0.329 | 0.522 | 0.730 | 0.405 | 2.527 | 2.511 | 2.538 |  |  |  | 0.96 |
| Sc5_d1 | CAHB | -0.437 | 0.419 | 0.546 | 0.910 | 0.165 | 2.203 |  | 2.203 |  | 0.47 | 102.7 | 0.81 |
| Sc5_d1 | CAHB-p10 | -0.290 | 0.581 | 0.549 | 0.975 | 0.170 | 2.428 |  | 2.428 |  | 0.32 | 66.5 | 0.13 |
| Sc5_d1 | CAHB-lit | -0.089 | 12.198 | 0.948 | 1.000 | 0.065 | 4.617 |  | 4.617 |  | 0.01 | 0.5 | 0.02 |
| Sc5_d1 | LRC-BART-100 | -0.025 | 0.191 | 0.210 | 0.900 | 1.000 | 0.772 |  | 0.772 |  | 0.10 | 0.9 | 2.14 |
| Sc5_d1 | BART-PP | -0.025 | 0.208 | 0.202 | 0.950 | 1.000 | 0.773 |  | 0.773 |  |  |  | 1.53 |
| Sc5_d1 | BART-CP | -0.700 | 0.147 | 0.716 | 0.000 | 0.515 | 1.038 |  | 1.038 |  |  |  | 1.59 |
| Sc5_d1 | BART-NP | -0.016 | 0.202 | 0.223 | 0.895 | 1.000 | 1.160 |  | 1.160 |  |  |  | 1.05 |
| Sc5_d1 | PSCL | -0.519 | 0.269 | 0.581 | 0.490 | 0.400 | 2.543 |  | 2.543 |  |  |  | 0.00 |
| Sc5_d1 | MAP-100 | -0.377 | 0.336 | 0.504 | 0.775 | 0.445 | 2.523 |  | 2.523 |  |  |  | 0.93 |
| Sc5_d2 | CAHB | -0.067 | 0.479 | 0.400 | 0.950 | 0.515 | 1.983 |  | 1.983 |  | 0.03 | 5.8 | 0.73 |
| Sc5_d2 | CAHB-p10 | -0.028 | 0.648 | 0.340 | 1.000 | 0.205 | 2.344 |  | 2.344 |  | 0.00 | 0.0 | 0.13 |
| Sc5_d2 | CAHB-lit | -0.090 | 12.165 | 0.947 | 1.000 | 0.065 | 4.609 |  | 4.609 |  | 0.00 | 0.4 | 0.02 |
| Sc5_d2 | LRC-BART-100 | -0.019 | 0.189 | 0.205 | 0.910 | 1.000 | 0.777 |  | 0.777 |  | 0.01 | 0.2 | 2.16 |
| Sc5_d2 | BART-PP | -0.034 | 0.209 | 0.206 | 0.960 | 1.000 | 0.773 |  | 0.773 |  |  |  | 1.54 |
| Sc5_d2 | BART-CP | -1.402 | 0.153 | 1.410 | 0.000 | 0.745 | 1.615 |  | 1.615 |  |  |  | 1.62 |
| Sc5_d2 | BART-NP | -0.016 | 0.202 | 0.223 | 0.895 | 1.000 | 1.160 |  | 1.160 |  |  |  | 1.05 |
| Sc5_d2 | PSCL | -1.021 | 0.269 | 1.054 | 0.040 | 0.050 | 2.689 |  | 2.689 |  |  |  | 0.00 |
| Sc5_d2 | MAP-100 | -0.430 | 0.345 | 0.551 | 0.735 | 0.370 | 2.533 |  | 2.533 |  |  |  | 0.94 |

### Table B. Paired differences from LRC-BART-100 over the same 200 replicates

rmse_diff is RMSE(method) minus RMSE(LRC-BART-100); its SE is the paired SE of the mean squared-error difference divided by the sum of the two RMSEs (delta method). The other columns are paired mean differences with their SEs (SD, coverage and power differences; surface-RMSE differences overall and by region).

| scenario | method | rmse | rmse_lrc | rmse_diff | rmse_diff_se | bias_diff | bias_diff_se | sd_diff | cov_diff | cov_diff_se | pow_diff | pow_diff_se | surf_diff | surfR_diff | surfRc_diff |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | CAHB | 0.258 | 0.175 | 0.082 | 0.015 | 0.006 | 0.014 | 0.240 | 0.065 | 0.017 | -0.220 | 0.029 | 1.419 | 1.331 | 1.501 |
| Sc1 | CAHB-p10 | 0.262 | 0.175 | 0.087 | 0.015 | -0.010 | 0.015 | 0.372 | 0.065 | 0.017 | -0.665 | 0.033 | 1.651 | 1.558 | 1.738 |
| Sc1 | CAHB-lit | 0.950 | 0.175 | 0.775 | 0.643 | -0.078 | 0.065 | 12.055 | 0.065 | 0.017 | -0.935 | 0.017 | 3.870 | 4.086 | 2.389 |
| Sc1 | BART-PP | 0.206 | 0.175 | 0.031 | 0.005 | -0.006 | 0.005 | 0.037 | 0.005 | 0.018 | 0.000 | 0.000 | 0.016 | 0.014 | 0.018 |
| Sc1 | BART-CP | 0.148 | 0.175 | -0.027 | 0.007 | 0.003 | 0.007 | -0.025 | 0.025 | 0.017 | 0.000 | 0.000 | -0.006 | -0.009 | -0.003 |
| Sc1 | BART-NP | 0.223 | 0.175 | 0.048 | 0.007 | -0.007 | 0.007 | 0.032 | -0.040 | 0.021 | 0.000 | 0.000 | 0.403 | 0.407 | 0.398 |
| Sc1 | PSCL | 0.262 | 0.175 | 0.086 | 0.014 | -0.009 | 0.014 | 0.098 | 0.005 | 0.022 | -0.025 | 0.011 | 1.734 | 1.616 | 1.846 |
| Sc1 | MAP-100 | 0.266 | 0.175 | 0.091 | 0.014 | -0.013 | 0.015 | 0.119 | 0.015 | 0.021 | -0.055 | 0.016 | 1.735 | 1.617 | 1.845 |
| Sc4_X5_d1 | CAHB | 0.333 | 0.227 | 0.105 | 0.018 | -0.143 | 0.017 | 0.219 | 0.110 | 0.024 | -0.575 | 0.035 | 1.353 | 1.356 | 1.351 |
| Sc4_X5_d1 | CAHB-p10 | 0.380 | 0.227 | 0.152 | 0.019 | -0.196 | 0.018 | 0.351 | 0.130 | 0.024 | -0.875 | 0.023 | 1.603 | 1.606 | 1.601 |
| Sc4_X5_d1 | CAHB-lit | 0.948 | 0.227 | 0.721 | 0.613 | -0.024 | 0.066 | 12.013 | 0.130 | 0.024 | -0.935 | 0.017 | 3.801 | 4.042 | 2.303 |
| Sc4_X5_d1 | BART-PP | 0.203 | 0.227 | -0.025 | 0.005 | 0.041 | 0.004 | 0.015 | 0.090 | 0.021 | 0.000 | 0.000 | -0.017 | -0.002 | -0.030 |
| Sc4_X5_d1 | BART-CP | 0.382 | 0.227 | 0.154 | 0.012 | -0.287 | 0.012 | -0.048 | -0.555 | 0.040 | -0.005 | 0.005 | 0.077 | -0.028 | 0.165 |
| Sc4_X5_d1 | BART-NP | 0.223 | 0.227 | -0.004 | 0.007 | 0.049 | 0.007 | 0.008 | 0.025 | 0.023 | 0.000 | 0.000 | 0.340 | 0.368 | 0.313 |
| Sc4_X5_d1 | PSCL | 0.373 | 0.227 | 0.146 | 0.017 | -0.204 | 0.016 | 0.076 | -0.005 | 0.030 | -0.235 | 0.030 | 1.685 | 1.671 | 1.701 |
| Sc4_X5_d1 | MAP-100 | 0.381 | 0.227 | 0.153 | 0.017 | -0.189 | 0.017 | 0.109 | 0.020 | 0.028 | -0.285 | 0.032 | 1.685 | 1.664 | 1.706 |
| Sc4_X5_d2 | CAHB | 0.453 | 0.206 | 0.247 | 0.022 | -0.313 | 0.018 | 0.234 | 0.050 | 0.025 | -0.750 | 0.031 | 1.323 | 1.352 | 1.295 |
| Sc4_X5_d2 | CAHB-p10 | 0.533 | 0.206 | 0.327 | 0.026 | -0.245 | 0.027 | 0.396 | 0.070 | 0.023 | -0.830 | 0.027 | 1.576 | 1.566 | 1.586 |
| Sc4_X5_d2 | CAHB-lit | 0.947 | 0.206 | 0.742 | 0.624 | -0.069 | 0.066 | 11.999 | 0.095 | 0.021 | -0.935 | 0.017 | 3.771 | 3.995 | 2.293 |
| Sc4_X5_d2 | BART-PP | 0.203 | 0.206 | -0.003 | 0.004 | -0.009 | 0.003 | 0.022 | 0.045 | 0.016 | 0.000 | 0.000 | -0.035 | -0.039 | -0.029 |
| Sc4_X5_d2 | BART-CP | 0.714 | 0.206 | 0.508 | 0.016 | -0.677 | 0.010 | -0.039 | -0.900 | 0.022 | -0.475 | 0.035 | 0.395 | -0.048 | 0.711 |
| Sc4_X5_d2 | BART-NP | 0.223 | 0.206 | 0.018 | 0.007 | 0.005 | 0.006 | 0.014 | -0.010 | 0.022 | 0.000 | 0.000 | 0.316 | 0.327 | 0.307 |
| Sc4_X5_d2 | PSCL | 0.581 | 0.206 | 0.375 | 0.024 | -0.499 | 0.015 | 0.084 | -0.410 | 0.042 | -0.615 | 0.034 | 1.699 | 1.743 | 1.658 |
| Sc4_X5_d2 | MAP-100 | 0.522 | 0.206 | 0.316 | 0.025 | -0.390 | 0.018 | 0.140 | -0.175 | 0.039 | -0.595 | 0.035 | 1.683 | 1.692 | 1.676 |
| Sc5_d1 | CAHB | 0.546 | 0.210 | 0.335 | 0.023 | -0.413 | 0.019 | 0.228 | 0.010 | 0.027 | -0.835 | 0.026 | 1.431 |  | 1.431 |
| Sc5_d1 | CAHB-p10 | 0.549 | 0.210 | 0.339 | 0.025 | -0.266 | 0.027 | 0.390 | 0.075 | 0.023 | -0.830 | 0.027 | 1.655 |  | 1.655 |
| Sc5_d1 | CAHB-lit | 0.948 | 0.210 | 0.738 | 0.622 | -0.064 | 0.066 | 12.007 | 0.100 | 0.021 | -0.935 | 0.017 | 3.845 |  | 3.845 |
| Sc5_d1 | BART-PP | 0.202 | 0.210 | -0.009 | 0.003 | -0.000 | 0.003 | 0.017 | 0.050 | 0.015 | 0.000 | 0.000 | 0.001 |  | 0.001 |
| Sc5_d1 | BART-CP | 0.716 | 0.210 | 0.506 | 0.016 | -0.676 | 0.010 | -0.044 | -0.900 | 0.021 | -0.485 | 0.035 | 0.266 |  | 0.266 |
| Sc5_d1 | BART-NP | 0.223 | 0.210 | 0.013 | 0.006 | 0.009 | 0.006 | 0.011 | -0.005 | 0.021 | 0.000 | 0.000 | 0.387 |  | 0.387 |
| Sc5_d1 | PSCL | 0.581 | 0.210 | 0.371 | 0.024 | -0.495 | 0.016 | 0.078 | -0.410 | 0.042 | -0.600 | 0.035 | 1.771 |  | 1.771 |
| Sc5_d1 | MAP-100 | 0.504 | 0.210 | 0.294 | 0.025 | -0.353 | 0.018 | 0.145 | -0.125 | 0.036 | -0.555 | 0.035 | 1.751 |  | 1.751 |
| Sc5_d2 | CAHB | 0.400 | 0.205 | 0.195 | 0.038 | -0.048 | 0.022 | 0.290 | 0.040 | 0.025 | -0.485 | 0.035 | 1.206 |  | 1.206 |
| Sc5_d2 | CAHB-p10 | 0.340 | 0.205 | 0.134 | 0.019 | -0.009 | 0.018 | 0.459 | 0.090 | 0.020 | -0.795 | 0.029 | 1.566 |  | 1.566 |
| Sc5_d2 | CAHB-lit | 0.947 | 0.205 | 0.742 | 0.624 | -0.071 | 0.066 | 11.977 | 0.090 | 0.020 | -0.935 | 0.017 | 3.832 |  | 3.832 |
| Sc5_d2 | BART-PP | 0.206 | 0.205 | 0.001 | 0.003 | -0.015 | 0.003 | 0.020 | 0.050 | 0.015 | 0.000 | 0.000 | -0.004 |  | -0.004 |
| Sc5_d2 | BART-CP | 1.410 | 0.205 | 1.205 | 0.019 | -1.383 | 0.010 | -0.036 | -0.910 | 0.020 | -0.255 | 0.031 | 0.838 |  | 0.838 |
| Sc5_d2 | BART-NP | 0.223 | 0.205 | 0.018 | 0.006 | 0.003 | 0.006 | 0.013 | -0.015 | 0.021 | 0.000 | 0.000 | 0.383 |  | 0.383 |
| Sc5_d2 | PSCL | 1.054 | 0.205 | 0.849 | 0.029 | -1.002 | 0.015 | 0.080 | -0.870 | 0.026 | -0.950 | 0.015 | 1.912 |  | 1.912 |
| Sc5_d2 | MAP-100 | 0.551 | 0.205 | 0.346 | 0.029 | -0.411 | 0.019 | 0.156 | -0.175 | 0.037 | -0.630 | 0.034 | 1.755 |  | 1.755 |


## 6. What CAHB does relative to LRC-BART, scenario by scenario

Sc1 (compatible RWD). CAHB borrows on every replicate (map 0.52 in both regions, realized ESS median 113, 10% to 90% range 93 to 140) and is unbiased (-0.003). Its RMSE is 0.258 against 0.175 for LRC-BART-100 (paired difference 0.082, SE 0.015), 0.206 for BART-PP and 0.223 for BART-NP; it sits with the marginal-mean methods PSCL (0.262) and MAP-100 (0.266). The posterior SD is 0.41, 1.6 times the MC RMSE, so coverage is 1.00 and power 0.78 (LRC-BART-100: 1.00). The control-surface RMSE is 2.18 against 0.76: the four-covariate kernel at the local-sample-20 bandwidth smooths the surface almost to its mean (the no-borrowing surface RMSE at this bandwidth is 2.2 to 2.5, and the sample SD of the true surface is about 2.5).

Sc4_X5_d1 and Sc4_X5_d2 (shift of 1 or 2 outside R). CAHB does not localize: the map is 0.53 in R and 0.48 in R^c at delta 1, and 0.49 against 0.41 at delta 2, where LRC-BART-100 has 0.64 against 0.53 and 0.66 against 0.05. Borrowing in R^c produces bias -0.21 at delta 1 and -0.33 at delta 2, RMSE 0.333 and 0.453 against 0.227 and 0.206 for LRC-BART-100 (paired differences 0.105, SE 0.018, and 0.247, SE 0.022); coverage 0.98 and 0.96 (the wide posterior absorbs the bias), power 0.43 and 0.25. At delta 2 the switch-off is all-or-nothing on 7.5% of replicates (ESS 0). The surface RMSE in R is 2.17 against 0.82. The bandwidth needed for 20 local subjects in four dimensions is about 1.4 units of X5, which is why the X5 = 2 boundary is not resolved; a one-dimensional kernel in X5 would resolve it, but that would require knowing which covariate carries the shift.

Sc5_d1 (global shift of 1). CAHB misses a shift of one residual SD on 94% of replicates (map 0.47, ESS median 108) and is biased by -0.44, RMSE 0.546 against 0.210 (paired difference 0.335, SE 0.023); coverage 0.91, power 0.17. On the 12 replicates where the switch fired the bias was +0.38. Among the reference methods only BART-CP (-0.70) and PSCL (-0.52) are more biased; BART-PP (-0.03), BART-NP (-0.02) and LRC-BART-100 (-0.03) are not.

Sc5_d2 (global shift of 2). The tau switch fires on 95% of replicates: map 0.03, ESS 0, bias -0.07, RMSE 0.400 against 0.205 (paired difference 0.195, SE 0.038); on the 10 replicates where it did not fire the bias was -1.15. The RMSE without borrowing (0.31 on the 190 switched-off replicates) is what the kernel estimator alone delivers, which is why CAHB's RMSE is above LRC-BART's even when it borrows nothing.

CAHB-p10 shares the pattern with a larger posterior SD (0.54 to 0.65), lower power (0.34 in Sc1) and a more frequent all-or-nothing switch (39% of Sc4_X5_d2 and 36% of Sc5_d1 replicates at ESS 0). CAHB-lit is degenerate: the interpolating kernel gives phi_0 about 0.07, a posterior SD of 12, ESS 0.4, coverage 1.00 and power 0.065 in every scenario, and a surface RMSE of 4.6.

In sum, on this design CAHB borrows about 110 effective controls when the RWD is compatible, switches off at a shift of two residual SDs and not at one, does not localize a regional shift, and carries the kernel smoother's own error (surface RMSE above 2, posterior SD 0.4) in every scenario. Relative to LRC-BART-100 its RMSE is higher by 0.08 (Sc1), 0.11 and 0.25 (Sc4 delta 1, 2), 0.34 and 0.20 (Sc5 delta 1, 2), with paired SEs of 0.015 to 0.038. These results are for the transfer of CAHB to ten continuous covariates with a nonlinear surface; the published setting has four covariates, two of them binary, and a near-linear surface.

## 7. Bib-ready description for the comparator list

Jin2023CAHB (already in `05-writing/merged.bib`). Suggested sentence: "CAHB \citep{Jin2023CAHB} places a local prior $\mu_0(x) \sim N\{\theta_0(x), 1/\tau(x)\}$ on the current control mean, with $\theta_0$ a historical mean function supplied as an input and a kernel-local precision $\tau(x)$ that is thresholded and $L_1$-projected so that borrowing is switched off where the local discrepancy between the sources is large; we transfer it to our design with a four-covariate kernel, a local sample of 20 and the published constants matched on the dimensionless scale (Web Appendix)."

## 8. Unresolved

The transfer choices in Section 2 (kernel covariates, bandwidth rule, the same-kernel theta_0, the scale matching of gamma and lambda_2, the joint draws) are ours; the paper offers no rule for any of them beyond its own four-covariate setting. The literal transfer is reported so the reader can see that it is degenerate rather than take our word for it. If a reviewer prefers the authors' independent per-profile draws for the primary method, the point estimates, surfaces and maps are unchanged and only the SD, coverage and power columns would change (`fit_cahb(..., joint_draws = FALSE)`). The survival study was not attempted.
