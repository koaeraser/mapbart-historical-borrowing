# Phase 2 development go/no-go run (Gaussian, R = 200)

Date 2026-09-09. Run: 12 scenarios, 16 methods, R = 200, seed 2026, 16 forked workers; 2400 replicate jobs, no fit errors; wall time 50.8 min (fits 50.5 min, from 11:37:58 to 12:28:44). Code: `01-code` after the post-verification fixes (`lrcbart/NOTES.md`, "Post-verification fixes"). Files here: `run_dev.R`, the `summarize_sim()` CSVs, `tables.md` (all twelve scenario tables), `criteria.csv`, `wall_time.txt`.

## Verdict

GO-WITH-CHANGES. The sampler and the calibration behave as designed and the method wins on efficiency (Sc1), robustness (Sc3, on RMSE and coverage) and abstention (Sc5), and it dominates C-BART and the MAP comparators everywhere they matter. It does not meet the Sc4 local-borrowing criteria: the control surface in the compatible region is not better than BART-PP's, the map in R is 0.65 against the restated 0.7, and the learned discrepancy under-shoots the true shift. Ten of twenty-five criteria fail. The changes are listed at the end.

## Criterion table (LRC-BART-100, spec Section E as written; items 16, 17 and 25 are the restated map criterion)

| # | scenario | criterion | value | result |
|---|---|---|---|---|
| 1 | Sc3_rho-0.5 | abs(bias) <= 0.05 | 0.064 | FAIL |
| 2 | Sc3_rho-0.5 | coverage >= 0.92 | 0.930 | PASS |
| 3 | Sc3_rho-0.5 | RMSE <= BART-PP + 0.02 | 0.270 vs 0.261 + 0.02 | PASS |
| 4 | Sc3_rho0 | abs(bias) <= 0.05 | 0.062 | FAIL |
| 5 | Sc3_rho0 | coverage >= 0.92 | 0.940 | PASS |
| 6 | Sc3_rho0 | RMSE <= BART-PP + 0.02 | 0.269 vs 0.262 + 0.02 | PASS |
| 7 | Sc3_rho0.5 | abs(bias) <= 0.05 | 0.063 | FAIL |
| 8 | Sc3_rho0.5 | coverage >= 0.92 | 0.955 | PASS |
| 9 | Sc3_rho0.5 | RMSE <= BART-PP + 0.02 | 0.266 vs 0.254 + 0.02 | PASS |
| 10 | Sc1 | RMSE <= 0.90 x BART-PP | 0.170 vs 0.185 | PASS |
| 11 | Sc1 | power >= BART-PP + 0.05 | 1.000 vs 1.000 + 0.05 | FAIL |
| 12 | Sc4_X5_d2 | abs(bias) <= 0.10 | -0.026 | PASS |
| 13 | Sc4_X5_d2 | control RMSE in R <= 0.90 x BART-PP's | 0.815 vs 0.702 | FAIL |
| 14 | Sc4_X5_d2 | all-spike map in R >= 0.5 | 0.296 | FAIL |
| 15 | Sc4_X5_d2 | all-spike map in R^c <= 0.15 | 0.023 | PASS |
| 16 | Sc4_X5_d2 | P(abs(g) < 0.5) in R >= 0.7 (restated) | 0.654 | FAIL |
| 17 | Sc4_X5_d2 | P(abs(g) < 0.5) in R^c <= 0.15 (restated) | 0.057 | PASS |
| 18 | Sc4_X5_d2 | mean g within 0.25 of 0 in R | -0.262 | FAIL |
| 19 | Sc4_X5_d2 | mean g within 0.25 of -2 in R^c | -1.682 | FAIL |
| 20 | Sc4_X5_d1 | abs(bias) <= 0.10 | -0.071 | PASS |
| 21 | Sc4_X5_d1 | overall RMSE no worse than BART-PP | 0.230 vs 0.203 | FAIL |
| 22 | Sc5_d2 | abs(bias) <= 0.10 | -0.027 | PASS |
| 23 | Sc5_d2 | RMSE <= 1.10 x BART-NP | 0.203 vs 0.246 | PASS |
| 24 | Sc5_d2 | all-spike map <= 0.15 | 0.004 | PASS |
| 25 | Sc5_d2 | P(abs(g) < 0.5) <= 0.15 (restated) | 0.009 | PASS |

Most likely reasons for the failures. Items 1, 4, 7: the bias of 0.06 (Monte Carlo SE 0.02; BART-PP 0.036, the free-offset LRC-BART-w0 0.030) comes from the residual spike mass under a global shift (map 0.08 to 0.12, realized ESS 3), which pulls the constant discrepancy of 1.3 slightly toward zero; the criterion is missed by 0.01 while coverage and RMSE pass. Item 11: the harness's power is rejection of zero by the 95% interval at delta = 1, which is 1.00 for every BART method, so the criterion cannot discriminate (the paper's power column is a 90% interval at an effect of 0.5; see `01-code/reports/reference_check.md`). Items 13, 16, 18, 19, 21: the g ensemble under-shoots the shift (mean g -1.68 in R^c against -2, and -0.26 in R against 0, that is the shift leaks across the X5 = 2 boundary), so part of the R^c discrepancy is absorbed in R and part of the shift is left in f; the paired surface RMSE difference in R is +0.035 (SE 0.006) against BART-PP, and the leak is also the -0.07 ATE bias at delta = 1. BART-PP, with the source indicator as a covariate in one 50-tree forest, absorbs a regional shift with less error than a separate 10-tree g ensemble with the (0.5, 3) depth prior and a fixed slab scale. Item 14: the all-spike map is an H_g-dependent conjunction, as the verifier found.

## Scenario tables (bias, SD, RMSE, coverage, power of the ATE; Sc4 columns add the control-surface RMSE in R and R^c, the new map, the all-spike map, the mean g, and the prior ESS_0, realized ESS, ESS ceiling and capped fraction). All twelve are in `tables.md`; the decisive ones follow.

### Sc1 (best linear method LM-CP)

| method | bias | SD | RMSE | cov | power | map | ESS0 | ESSreal | ceiling | capped |
|---|---|---|---|---|---|---|---|---|---|---|
| LM-CP | -0.003 | 0.243 | 0.236 | 0.95 | 0.99 | | | | | |
| MAP-100 | -0.022 | 0.290 | 0.266 | 0.95 | | | 100 | | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | 1.00 | | | | | |
| BART-CP | -0.006 | 0.145 | 0.148 | 0.96 | 1.00 | | | | | |
| BART-PP | -0.014 | 0.207 | 0.206 | 0.94 | 1.00 | | | | | |
| C-BART | -0.007 | 0.154 | 0.151 | 0.96 | 1.00 | 1.00 | 81.9 | 81.3 | 121 | 0.23 |
| LRC-BART-w0 | -0.011 | 0.189 | 0.209 | 0.91 | 1.00 | 0.69 | 81.6 | 0.1 | 123 | 0.22 |
| LRC-BART-50 | -0.007 | 0.176 | 0.181 | 0.94 | 1.00 | 0.92 | 42.6 | 30.4 | 123 | 0.01 |
| LRC-BART-75 | -0.008 | 0.172 | 0.178 | 0.93 | 1.00 | 0.92 | 62.1 | 42.4 | 123 | 0.06 |
| LRC-BART-100 | -0.009 | 0.170 | 0.170 | 0.94 | 1.00 | 0.93 | 81.7 | 55.2 | 123 | 0.24 |

### Sc3, rho = 0 (best linear method LM-PP; rho = -0.5 and 0.5 are within 0.01 of these)

| method | bias | SD | RMSE | cov | power | map | ESS0 | ESSreal |
|---|---|---|---|---|---|---|---|---|
| LM-PP | 0.051 | 0.356 | 0.332 | 0.96 | 0.86 | | | |
| MAP-100 | 0.372 | 0.358 | 0.501 | 0.81 | | | 100 | |
| BART-NP | 0.008 | 0.291 | 0.276 | 0.95 | 0.94 | | | |
| BART-CP | 0.980 | 0.202 | 1.000 | 0.00 | 1.00 | | | |
| BART-PP | 0.037 | 0.276 | 0.262 | 0.95 | 0.97 | | | |
| C-BART | 0.237 | 0.299 | 0.380 | 0.90 | 0.98 | 0.15 | 84.8 | 22.5 |
| LRC-BART-w0 | 0.030 | 0.279 | 0.262 | 0.96 | 0.97 | 0.10 | 84.3 | 0.2 |
| LRC-BART-50 | 0.066 | 0.282 | 0.273 | 0.95 | 0.97 | 0.08 | 44.3 | 3.1 |
| LRC-BART-100 | 0.062 | 0.282 | 0.269 | 0.94 | 0.97 | 0.08 | 84.1 | 2.7 |

### Sc4, region X5 > 2, delta_rwd = 2 (best linear method LM-NP)

| method | bias | SD | RMSE | cov | RMSE_R | RMSE_Rc | map_R | map_Rc | spike_R | spike_Rc | g_R | g_Rc | ESS0 | ESSreal | ceiling | capped |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LM-NP | -0.020 | 0.340 | 0.339 | 0.93 | 2.433 | 2.339 | | | | | | | | | | |
| MAP-100 | -0.410 | 0.329 | 0.522 | 0.73 | | | | | | | | | 100 | | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | 1.146 | 1.169 | | | | | | | | | | |
| BART-CP | -0.697 | 0.149 | 0.714 | 0.01 | 0.772 | 1.573 | | | | | | | | | | |
| BART-PP | -0.029 | 0.210 | 0.203 | 0.95 | 0.780 | 0.832 | | | | | | | | | | |
| C-BART | -0.247 | 0.200 | 0.408 | 0.66 | 0.801 | 1.116 | 0.63 | 0.41 | 1.00 | 1.00 | -0.34 | -0.91 | 81.7 | 45.6 | 120 | 0.24 |
| LRC-BART-w0 | -0.020 | 0.191 | 0.206 | 0.92 | 0.831 | 0.874 | 0.55 | 0.05 | 0.00 | 0.00 | -0.21 | -1.75 | 81.4 | 0.1 | 121 | 0.24 |
| LRC-BART-50 | -0.026 | 0.189 | 0.208 | 0.90 | 0.816 | 0.869 | 0.66 | 0.06 | 0.31 | 0.04 | -0.27 | -1.67 | 42.4 | 1.5 | 119 | 0.01 |
| LRC-BART-75 | -0.021 | 0.188 | 0.203 | 0.90 | 0.824 | 0.865 | 0.64 | 0.05 | 0.29 | 0.02 | -0.27 | -1.68 | 61.2 | 1.4 | 121 | 0.05 |
| LRC-BART-100 | -0.026 | 0.188 | 0.206 | 0.90 | 0.815 | 0.863 | 0.65 | 0.06 | 0.30 | 0.02 | -0.26 | -1.68 | 81.7 | 1.5 | 120 | 0.26 |

Sc4 with delta = 1 (X5): LRC-BART-100 bias -0.071, RMSE 0.230, coverage 0.88, RMSE_R 0.774 (BART-PP 0.776), map 0.67 / 0.55, g -0.29 / -0.50. The X7 region reproduces the X5 numbers within 0.01 at both deltas.

### Sc5, delta_rwd = 2 (best linear method LM-NP)

| method | bias | SD | RMSE | cov | power | map | ESSreal |
|---|---|---|---|---|---|---|---|
| LM-NP | -0.020 | 0.340 | 0.339 | 0.93 | 0.81 | | |
| MAP-100 | -0.430 | 0.345 | 0.551 | 0.74 | | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | 1.00 | | |
| BART-CP | -1.402 | 0.153 | 1.410 | 0.00 | 0.74 | | |
| BART-PP | -0.034 | 0.209 | 0.206 | 0.96 | 1.00 | | |
| C-BART | -0.318 | 0.198 | 0.643 | 0.74 | 0.94 | 0.18 | 43.7 |
| LRC-BART-w0 | -0.016 | 0.189 | 0.208 | 0.92 | 1.00 | 0.01 | 0.1 |
| LRC-BART-100 | -0.027 | 0.189 | 0.203 | 0.93 | 1.00 | 0.01 | 0.5 |

## Key comparisons

Efficiency (Sc1). LRC-BART-100 has RMSE 0.170 against BART-PP's 0.206 (paired squared-error difference -0.013, SE 0.002), with coverage 0.94; C-BART's 0.151 is the floor and BART-CP's 0.148. The prior ESS_0 is 82 at the requested 100 because the Pr-rule chooses the s_0^2 at which 95 percent of blocks are at or below the target; the realized ESS is 55.

Robustness (Sc3). RMSE 0.266 to 0.270 against BART-PP's 0.254 to 0.262 and BART-NP's 0.266 to 0.276, coverage 0.93 to 0.955, bias 0.06 (BART-PP 0.036). C-BART carries bias 0.24 to 0.30, MAP-100 0.37, BART-CP 0.84 to 0.98.

Local borrowing (Sc4). The map separates the regions (0.65 in R, 0.06 in R^c at delta 2, on both X5 and X7), the ATE bias is -0.03 and coverage 0.90, and the model abstains in R^c (spike map 0.02). It does not gain in R: the surface RMSE in R is 0.815 against BART-PP's 0.780 and LRC-BART-w0's 0.831, so the spike improves on the free offset by 0.016 (SE 0.006) but the f + g parameterization costs 0.05 against the single forest. At delta 1 the ATE RMSE is 0.230 against 0.203.

Abstention (Sc5). RMSE 0.203 against BART-NP's 0.223 and BART-PP's 0.206, bias -0.03, map 0.01, realized ESS 0.5: the method abstains and loses nothing to BART-NP.

Does the slab matter (C-BART vs LRC-BART). Yes: Sc4 delta 2 RMSE 0.408 vs 0.206 (bias -0.25 vs -0.03), Sc5 delta 2 0.643 vs 0.203, Sc3 0.38 to 0.44 vs 0.27; in Sc1 the two are within 0.02.

Sensitivity to N_target. Sc1 RMSE 0.181 / 0.178 / 0.170 for 50 / 75 / 100; elsewhere the three variants are within 0.005 of one another. In Sc2 all three are capped in every replicate (ceiling 22 in RCT units, prior ESS_0 24) and give the same model, RMSE 0.185 to 0.188 against BART-PP's 0.198. In Sc1 the ceiling is 73 to 176 across replicates (median 122) and N_target = 100 is capped in 24 percent.

MAP comparator. MAP-100 is worse than LRC-BART-100 in every scenario: Sc1 0.266 vs 0.170, Sc3 0.501 vs 0.269, Sc4 delta 2 0.522 vs 0.206, Sc5 0.551 vs 0.203; its bias in Sc3 to Sc5 is 0.37 to 0.43 with coverage 0.73 to 0.81.

## Changes before the full run

1. Sc4 mechanism. The discrepancy ensemble under-shoots the shift and leaks across the region boundary. Test, on Sc4 at delta 2 with R about 50: a less restrictive g tree prior (alpha_g 0.95, beta_g 2, or n_min 2 for g), a larger slab scale tau_1^2, and H_g = 5 with the same changes; target g_Rc within 0.25 of -2 and g_R within 0.25 of 0, then re-check Sc1 and Sc3.
2. Surface criterion. Replace "RMSE in R at most 0.90 times BART-PP's" by a paired comparison with its Monte Carlo SE; the R-region surface error is dominated by f's approximation error common to both methods, and a 10 percent margin is not reachable by borrowing alone.
3. ESS ladder. Report the ESS ceiling sigma_1^2 / V_mu^f with every fit and define the targets relative to it (for instance 25, 50 and 75 percent of the ceiling), since N_target = 100 is infeasible in 24 percent of Sc1 replicates and always in Sc2; keep the capped fraction as a column.
4. Map criterion. Adopt P(abs(g(x)) < c | data) with c pre-specified and drop the all-spike map from the criteria (keep it as a diagnostic); at c = 0.5 the R-region mean is 0.65, so either c or the threshold must be argued before the full run.
5. Power. Use the paper's definition (90 percent interval, effect 0.5) so that the Sc1 power criterion can discriminate.
6. Sc3 bias. Keep the 0.05 criterion and re-check after change 1; if the 0.06 persists, the manuscript states the bias against BART-PP's 0.036 rather than as an unqualified pass.
