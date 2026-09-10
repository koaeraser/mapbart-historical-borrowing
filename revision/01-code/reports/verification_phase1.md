# Phase 1 verification of the lrcbart sampler (independent, adversarial)

Date: 2026-09-09. The verifier did not write the code and modified nothing in the package or harness; check scripts and raw outputs are in the session scratchpad (`check1_leaf.R` to `check9_sc2.R`).

## Table of checks

| # | Check | Result | Number |
|---|---|---|---|
| 1 | M_f, M_g vs independent R closed form and brute-force integration (300 random leaves each, one-source, empty, n_1 = 0, extreme ybar) | PASS | max abs error 7e-14 both; n_1 = 0 gives log M_g = 0; P(z) vs T1 2e-16 |
| 2a | Grow/prune/change ratios against Chipman et al. 2010 / bartMachine | PASS | correct at `lrcbart.cpp` 323-328, 362-364, 375-377, 393 |
| 2b | Prior-predictive: single-arm g (no data, n_min = 0), 20 trees x 4000 draws vs forward simulation of the (0.5, 3) prior | PASS | P(stump) 0.4998 vs 0.5 (SE 0.0017); leaf-count distribution (1..5) 0.500/0.442/0.055/0.004/0.0001 vs 0.501/0.439/0.056/0.004/0.0002 |
| 3 | tau_0^2, tau_1^2, w, sigma_g^2 updates vs exact (truncated) densities by quadrature | PASS | tau_0^2 mean 0.009812 vs 0.009820, CDF at quantiles 0.0495/0.4989/0.9495; binding truncation 0.008367 vs 0.008366; w 0.80013 vs 0.80000; sigma_1^2 1.4065 vs 1.4077, sigma_2^2 3.7871 vs 3.7869 |
| 4 | Two-leaf exactness: H_f = 0, H_g = 1, binary covariate, (tree, z, theta, tau_0^2, w, sigma_1^2) all free, vs a 600 x 600 grid with (tree, z) enumerated; 400k draws | PASS | max abs error: posterior means of g 0.0007 (MCMC SE 0.0004), P(z) 0.0009 (SE 0.0005), P(split) 0.5377 vs 0.5379. A second configuration with P(split) = 0.9999: 0.0014 / 0.0021 |
| 5 | Plain BART vs `BART::wbart`, harness Sc1, 5 seeds, control-surface RMSE at the RCT profiles | PASS | pooled H_g = 0: 0.742 vs 0.723 (corr 0.988; the package keeps two residual variances); single-arm w = 1 on RWD: 0.784 vs 0.827; RCT-only: 1.031 vs 1.027 |
| 6a | Sc4 (delta = 2, R = {X5 > 2}), H_f = 50, H_g = 10, harness seeds 1-5 | CONCERN | ATE bias -0.004 (SD 0.18); control bias R +0.03, R^c +0.02; B_R 0.33 (0.00, 0.60, 0.64, 0.33, 0.07), B_Rc 0.03; per-leaf P(spike) 0.79, spike fraction 0.79; mean g R -0.18, R^c -1.87; RMSE_R 0.790 vs BART-PP 0.769 |
| 6b | Same at H_g = 5 | CONCERN | ATE bias -0.006; B_R 0.36, B_Rc 0.02; per-leaf P(spike) 0.70; g R -0.29, R^c -1.77; RMSE_R 0.788 |
| 6c | Alternative map P(abs(g(x)) < c) | reported | H_g = 10: c = 0.5 gives R 0.74, R^c 0.04; c = 1 gives R 0.90, R^c 0.06. H_g = 5: 0.73/0.06 and 0.87/0.08 |
| 7 | `time_methods` through `fit_one`, Sc2, 3 reps | PASS | Gaussian: LRC-BART-100 1.34 s, C-BART 1.30 s, BART-PP 0.94 s; survival: 1.52 s, 1.49 s, 1.25 s. NOTES' 1.3 s / 1.6 s hold |
| 8a | Wrapper defaults vs spec (H_g, alpha_g, beta_g, nu_0, a_w, b_w, tau_1^2 rule, n_min, k, nu_sigma, q) | PASS | all match |
| 8b | ESS calibration vs Section C | CONCERN | see bugs 1 to 3 |
| 8c | Package test suite | as documented | 1 failure of 143 expectations (test 7, B_R 0.451) |

## Bugs and concerns found

1. `R/lrcbart.R` lines 284-285 and `R/fit_lrcbart.R` line 42: the calibration sets sigma_1^2 to the Stage-1 RWD residual variance (3.0 to 3.2 in Sc1) although RCT controls are in hand; Section C allows this only "before RCT outcomes exist". The prior calibrated for N_target = 100 has ESS_0 = 58 to 62 with sigma_1^2 = 2.25 (three seeds), and the ceiling of 210 becomes 130 to 156. The realized ESS (`lrcbart.R` line 355) uses the posterior RCT sigma_1^2, so prior and realized ESS are in different units. Fix: pass an RCT residual estimate as `sigma1_sq`.
2. `R/lrcbart.R` lines 282-283: V_mu^f is the variance within blocks of 50 consecutive draws (lag-1 autocorrelation 0.2 to 0.3), 5 to 16 percent below the whole-chain variance (0.0156 vs 0.0186, 0.0156 vs 0.0171, 0.0180 vs 0.0189), inflating the ceiling and s_0^2 in the same direction as bug 1. Line 281 is dead code.
3. Feasibility under covariate shift. In Sc2 the ceiling sigma_1^2 / V_mu^f is 38 to 39 (V_mu^f 0.07 to 0.08 against 0.016 in Sc1), so N_target = 100 and 75 are infeasible: the code warns on every Sc2 replicate and returns s_0^2 = 1e-6, making LRC-BART-100 and LRC-BART-75 the same model (a point-mass spike) and LRC-BART-50 nearly so. This follows T3, but the N_target ladder collapses in Sc2; `at_boundary` is recorded in `ess`.
4. Sc4 go/no-go on 5 seeds: two of five criteria fail. B_R = 0.33 (criterion 0.5) confirms NOTES; the control-surface RMSE in R is 0.790 against 0.769 for BART-PP (criterion 0.69; paired difference +0.02, SE 0.04), which test 7 does not compute. BART-PP with the D indicator absorbs a regional shift as well as LRC-BART does, and the surface error is dominated by f's approximation error. These are method-level results, not sampler defects.
5. The map B(x) = P(all H_g leaves spike) is dominated by whole-ensemble slab draws: the fraction of draws in which B(x) = 0 for every RCT control x is 0.30 to 1.00 across seeds (seed 1: 0.999, B_R = 0.00). H_g = 5 does not repair it (B_R 0.36).
6. Spec text: "K_0 of the order of 20 to 30" is inconsistent with the (0.5, 3) prior (1.56 leaves per tree, L_g = 14 at H_g = 10, observed 13.9), and "0.95^H_g = 0.6" assumes a per-leaf spike probability the data do not give (0.79 to 0.91). Minor: `fit_lrcbart.R` line 39, `seed + 1L` overflows at `.Machine$integer.max` (probability 2^-31).

## Which map definition is defensible

B(x) is a conjunction over H_g trees: its scale depends on H_g (w^{H_g} a priori) and one slab stump anywhere sets it to zero for every x, so it measures whether the whole ensemble is in the spike, not whether the surface is commensurate at x, and no fixed threshold applies without tuning H_g. P(abs(g(x)) < c | data) is a statement about the outcome-scale discrepancy at x, invariant to H_g, available from the draws the package already returns, and it separates the regions on all five seeds. We consider it the defensible map provided c is pre-specified on the outcome scale (one third of sigma_RCT, c = 0.5 here, or a minimal clinically relevant shift) and described as a compatibility probability; B(x) can stay as a diagnostic. T2 should then be restated for g(x), and the Sc4 criteria revised before Phase 2.

## Verdict

GATE PASSED for the sampler: leaf marginals, tree-move ratios, conjugate updates, the two-leaf joint posterior and the plain-BART limits are exact to Monte Carlo error. Before Phase 2, three decisions: (1) pass the RCT residual variance to `lrc_ess_calibrate` and use the whole-chain V_mu^f (bugs 1, 2), which changes every calibrated s_0^2; (2) accept or redesign the N_target ladder for Sc2 (bug 3); (3) choose the map definition and reset the Sc4 criteria (bugs 4, 5), which fail on the sampler as verified.
