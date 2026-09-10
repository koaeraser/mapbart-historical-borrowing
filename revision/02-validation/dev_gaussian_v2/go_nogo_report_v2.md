# Phase 2 development go/no-go run v2 (Gaussian, R = 200, after Phase 2b tuning)

Date 2026-09-09. Run: 12 scenarios, 18 methods, R = 200, seed 2026 (the v1 seed, so replicate data are paired with v1), 16 forked workers; 2400 replicate jobs, no fit errors; wall time 58.7 min (fits 58.5 min, 13:19:22 to 14:17:50). Configuration (spec A3, `lrcbart/NOTES.md` "Phase 2b tuning"): H_g = 5, split prior (0.5, 3), range-rule slab, n_min_g = 5, one g sweep, w ~ Beta(1, 1). Files: `run_dev_v2.R`, the `summarize_sim()` CSVs, `tables.md` (all twelve tables), `criteria.csv`, `v1_vs_v2.md`, `paired.csv`, `wall_time.txt`.

Two planned methods are absent: `run_dev_v2.R` runs LRC-BART-f90 and f50 but not f25, and not LRC-BART-Hg10 (the v1 configuration), so there is no within-run control at R = 200. The Hg10 comparison uses the v1 run (the same configuration on the same replicate data) and the Phase 2b confirmation stage (`tuning_phase2b/summary_confirm.csv`, R = 60). Both should be added to the full run.

## Verdict

GO-WITH-CHANGES, 18 of 25 criteria (v1: 15). The tuned configuration removes the Sc3 bias failure (0.049 to 0.050 against 0.062 to 0.064, a paired change of -0.013 with SE 0.002) and holds Sc1 efficiency and Sc5 abstention; the Sc4 local-borrowing criteria still fail, by the same margins as in v1, and NOTES.md shows why no g configuration meets them. The method is not better than BART-PP anywhere in Sc4: at delta 2 the two are tied on the ATE (paired squared-error difference 0.001, SE 0.001) and LRC-BART loses 0.039 (SE 0.006) on the control surface in R; at delta 1 it loses on the ATE (0.227 against 0.203, paired difference 0.011, SE 0.002) with coverage 0.87. The full study should proceed with this configuration, with Hg10 and f25 added, the paper's power definition, and a manuscript that states the Sc4 result as a boundary-local limitation rather than as local borrowing demonstrated against BART-PP.

## Criterion table (LRC-BART-100, spec Section E as written; 16, 17 and 25 are the restated map criterion)

| # | scenario | criterion | value | result |
|---|---|---|---|---|
| 1 | Sc3_rho-0.5 | abs(bias) <= 0.05 | 0.0503 (SE 0.019) | FAIL |
| 2 | Sc3_rho-0.5 | coverage >= 0.92 | 0.935 | PASS |
| 3 | Sc3_rho-0.5 | RMSE <= BART-PP + 0.02 | 0.268 vs 0.261 + 0.02 | PASS |
| 4 | Sc3_rho0 | abs(bias) <= 0.05 | 0.049 | PASS |
| 5 | Sc3_rho0 | coverage >= 0.92 | 0.960 | PASS |
| 6 | Sc3_rho0 | RMSE <= BART-PP + 0.02 | 0.268 vs 0.262 + 0.02 | PASS |
| 7 | Sc3_rho0.5 | abs(bias) <= 0.05 | 0.049 | PASS |
| 8 | Sc3_rho0.5 | coverage >= 0.92 | 0.955 | PASS |
| 9 | Sc3_rho0.5 | RMSE <= BART-PP + 0.02 | 0.260 vs 0.254 + 0.02 | PASS |
| 10 | Sc1 | RMSE <= 0.90 x BART-PP | 0.175 vs 0.185 | PASS |
| 11 | Sc1 | power >= BART-PP + 0.05 | 1.000 vs 1.000 + 0.05 | FAIL |
| 12 | Sc4_X5_d2 | abs(bias) <= 0.10 | -0.020 | PASS |
| 13 | Sc4_X5_d2 | control RMSE in R <= 0.90 x BART-PP's | 0.820 vs 0.702 | FAIL |
| 14 | Sc4_X5_d2 | all-spike map in R >= 0.5 | 0.279 | FAIL |
| 15 | Sc4_X5_d2 | all-spike map in R^c <= 0.15 | 0.022 | PASS |
| 16 | Sc4_X5_d2 | P(abs(g) < 0.5) in R >= 0.7 (restated) | 0.660 | FAIL |
| 17 | Sc4_X5_d2 | P(abs(g) < 0.5) in R^c <= 0.15 (restated) | 0.053 | PASS |
| 18 | Sc4_X5_d2 | mean g within 0.25 of 0 in R | -0.244 | PASS |
| 19 | Sc4_X5_d2 | mean g within 0.25 of -2 in R^c | -1.710 | FAIL |
| 20 | Sc4_X5_d1 | abs(bias) <= 0.10 | -0.064 | PASS |
| 21 | Sc4_X5_d1 | overall RMSE no worse than BART-PP | 0.227 vs 0.203 | FAIL |
| 22 | Sc5_d2 | abs(bias) <= 0.10 | -0.019 | PASS |
| 23 | Sc5_d2 | RMSE <= 1.10 x BART-NP | 0.205 vs 0.246 | PASS |
| 24 | Sc5_d2 | all-spike map <= 0.15 | 0.003 | PASS |
| 25 | Sc5_d2 | P(abs(g) < 0.5) <= 0.15 (restated) | 0.008 | PASS |

Item 1 misses by 0.0003 with a Monte Carlo SE of 0.019; the three Sc3 biases are one estimate (0.049 to 0.050) at the criterion's edge, so the manuscript states the bias against BART-PP's 0.036. Item 11 cannot discriminate (every BART method has power 1 at delta = 1 by the harness definition). Item 18 now passes (-0.244 against -0.262 in v1), but the change is 0.018 with SE 0.010, not a real gain. Items 13, 14, 16, 19, 21 are the Sc4 failures carried over from v1.

## Scenario tables (ATE bias, SD, RMSE, coverage; Sc4 adds the control-surface RMSE in R and R^c, the P(abs(g) < 0.5) map means, mean g, and the prior ESS_0, realized ESS, ceiling and capped fraction)

### Sc1 (best linear method LM-CP)

| method | bias | SD | RMSE | cov | map | ESS0 | ESSreal | ceiling | capped |
|---|---|---|---|---|---|---|---|---|---|
| LM-CP | -0.003 | 0.243 | 0.236 | 0.95 | | | | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | | | | | |
| BART-CP | -0.006 | 0.145 | 0.148 | 0.96 | | | | | |
| BART-PP | -0.014 | 0.207 | 0.206 | 0.94 | | | | | |
| C-BART | -0.008 | 0.154 | 0.153 | 0.94 | 1.00 | 81.9 | 81.0 | 121 | 0.23 |
| LRC-BART-w0 | -0.010 | 0.189 | 0.206 | 0.92 | 0.80 | 81.6 | 0.1 | 123 | 0.22 |
| LRC-BART-50 | -0.009 | 0.176 | 0.182 | 0.94 | 0.92 | 42.6 | 29.3 | 123 | 0.01 |
| LRC-BART-75 | -0.008 | 0.173 | 0.179 | 0.94 | 0.93 | 62.0 | 41.3 | 123 | 0.06 |
| LRC-BART-100 | -0.009 | 0.171 | 0.175 | 0.94 | 0.92 | 81.8 | 52.3 | 123 | 0.24 |
| LRC-BART-f50 | -0.006 | 0.174 | 0.179 | 0.93 | 0.92 | 51.3 | 34.1 | 122 | 0.00 |
| LRC-BART-f90 | -0.007 | 0.170 | 0.172 | 0.94 | 0.93 | 87.0 | 55.1 | 122 | 0.00 |

### Sc3, rho = 0 (best linear method LM-PP; rho = -0.5 and 0.5 within 0.01)

| method | bias | SD | RMSE | cov | map | ESS0 | ESSreal |
|---|---|---|---|---|---|---|---|
| LM-PP | 0.051 | 0.356 | 0.332 | 0.96 | | | |
| BART-NP | 0.008 | 0.291 | 0.276 | 0.95 | | | |
| BART-CP | 0.980 | 0.202 | 1.000 | 0.00 | | | |
| BART-PP | 0.037 | 0.276 | 0.262 | 0.95 | | | |
| C-BART | 0.235 | 0.301 | 0.375 | 0.89 | 0.15 | 84.7 | 22.2 |
| LRC-BART-w0 | 0.030 | 0.279 | 0.263 | 0.96 | 0.06 | 84.4 | 0.2 |
| LRC-BART-50 | 0.054 | 0.282 | 0.272 | 0.95 | 0.07 | 44.2 | 2.0 |
| LRC-BART-100 | 0.049 | 0.283 | 0.268 | 0.96 | 0.07 | 84.1 | 1.7 |
| LRC-BART-f90 | 0.044 | 0.281 | 0.266 | 0.95 | 0.07 | 159.6 | 1.3 |

### Sc4, region X5 > 2, delta_rwd = 2 (best linear method LM-NP)

| method | bias | SD | RMSE | cov | RMSE_R | RMSE_Rc | map_R | map_Rc | g_R | g_Rc | ESS0 | ESSreal | ceiling | capped |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LM-NP | -0.020 | 0.340 | 0.339 | 0.93 | 2.433 | 2.339 | | | | | | | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | 1.146 | 1.169 | | | | | | | | |
| BART-CP | -0.697 | 0.149 | 0.714 | 0.01 | 0.772 | 1.573 | | | | | | | | |
| BART-PP | -0.029 | 0.210 | 0.203 | 0.95 | 0.780 | 0.832 | | | | | | | | |
| C-BART | -0.248 | 0.196 | 0.412 | 0.64 | 0.793 | 1.103 | 0.68 | 0.40 | -0.29 | -0.96 | 81.7 | 45.8 | 120 | 0.24 |
| LRC-BART-w0 | -0.021 | 0.192 | 0.209 | 0.92 | 0.832 | 0.875 | 0.58 | 0.06 | -0.25 | -1.72 | 81.4 | 0.1 | 121 | 0.24 |
| LRC-BART-50 | -0.025 | 0.189 | 0.204 | 0.92 | 0.819 | 0.869 | 0.64 | 0.05 | -0.25 | -1.70 | 42.4 | 0.9 | 119 | 0.01 |
| LRC-BART-75 | -0.023 | 0.190 | 0.203 | 0.93 | 0.826 | 0.865 | 0.64 | 0.05 | -0.25 | -1.71 | 61.3 | 0.8 | 121 | 0.05 |
| LRC-BART-100 | -0.020 | 0.188 | 0.206 | 0.91 | 0.820 | 0.862 | 0.66 | 0.05 | -0.24 | -1.71 | 81.7 | 0.9 | 120 | 0.26 |
| LRC-BART-f50 | -0.026 | 0.189 | 0.204 | 0.92 | 0.819 | 0.867 | 0.66 | 0.06 | -0.25 | -1.70 | 50.1 | 0.9 | 120 | 0.00 |
| LRC-BART-f90 | -0.023 | 0.188 | 0.206 | 0.91 | 0.818 | 0.866 | 0.66 | 0.05 | -0.25 | -1.71 | 85.6 | 1.0 | 121 | 0.00 |

All-spike map for LRC-BART-100: 0.28 in R, 0.02 in R^c. The X7 region reproduces the X5 numbers within 0.01 at both deltas.

### Sc4, X5, delta_rwd = 1

| method | bias | RMSE | cov | RMSE_R | RMSE_Rc | map_R | map_Rc | g_R | g_Rc | ESS0 | ESSreal | capped |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BART-PP | -0.023 | 0.203 | 0.96 | 0.776 | 0.826 | | | | | | | |
| C-BART | -0.224 | 0.287 | 0.69 | 0.749 | 0.926 | 0.92 | 0.90 | -0.14 | -0.16 | 81.5 | 73.5 | 0.24 |
| LRC-BART-w0 | -0.017 | 0.203 | 0.93 | 0.790 | 0.835 | 0.53 | 0.41 | -0.34 | -0.63 | 82.1 | 0.1 | 0.23 |
| LRC-BART-100 | -0.064 | 0.227 | 0.87 | 0.778 | 0.856 | 0.64 | 0.53 | -0.30 | -0.52 | 81.2 | 21.6 | 0.23 |
| LRC-BART-f50 | -0.055 | 0.219 | 0.88 | 0.777 | 0.853 | 0.65 | 0.55 | -0.32 | -0.52 | 50.9 | 14.7 | 0.00 |

### Sc5, delta_rwd = 2 (best linear method LM-NP)

| method | bias | SD | RMSE | cov | map | ESSreal |
|---|---|---|---|---|---|---|
| LM-NP | -0.020 | 0.340 | 0.339 | 0.93 | | |
| BART-NP | -0.016 | 0.202 | 0.223 | 0.90 | | |
| BART-CP | -1.402 | 0.153 | 1.410 | 0.00 | | |
| BART-PP | -0.034 | 0.209 | 0.206 | 0.96 | | |
| C-BART | -0.350 | 0.188 | 0.693 | 0.73 | 0.21 | 48.2 |
| LRC-BART-w0 | -0.013 | 0.188 | 0.204 | 0.92 | 0.01 | 0.1 |
| LRC-BART-100 | -0.019 | 0.189 | 0.205 | 0.91 | 0.01 | 0.2 |

## v1 (H_g = 10, Beta(4,1)) against v2 (H_g = 5, Beta(1,1)), LRC-BART-100, same replicate data

| scenario | bias v1 / v2 | RMSE v1 / v2 | cov v1 / v2 | RMSE_R v1 / v2 | map_R v1 / v2 | map_Rc v1 / v2 | g_R v1 / v2 | g_Rc v1 / v2 | ESSreal v1 / v2 |
|---|---|---|---|---|---|---|---|---|---|
| Sc1 | -0.009 / -0.009 | 0.170 / 0.175 | 0.94 / 0.94 | 0.731 / 0.740 | 0.93 / 0.92 | | | | 55.2 / 52.3 |
| Sc2 | 0.029 / 0.029 | 0.185 / 0.185 | 0.94 / 0.95 | 0.754 / 0.753 | 0.88 / 0.87 | | | | 14.6 / 13.7 |
| Sc3_rho-0.5 | 0.064 / 0.050 | 0.270 / 0.268 | 0.93 / 0.94 | | 0.12 / 0.11 | | 1.20 / 1.22 | | 3.1 / 2.1 |
| Sc3_rho0 | 0.062 / 0.049 | 0.269 / 0.268 | 0.94 / 0.96 | | 0.08 / 0.07 | | 1.30 / 1.31 | | 2.7 / 1.7 |
| Sc3_rho0.5 | 0.063 / 0.049 | 0.266 / 0.260 | 0.95 / 0.95 | | 0.10 / 0.09 | | 1.22 / 1.24 | | 2.9 / 2.0 |
| Sc4_X5_d0.5 | -0.064 / -0.067 | 0.206 / 0.208 | 0.89 / 0.88 | 0.742 / 0.745 | 0.85 / 0.85 | 0.84 / 0.84 | -0.14 / -0.14 | -0.17 / -0.17 | 46.4 / 46.3 |
| Sc4_X5_d1 | -0.071 / -0.064 | 0.230 / 0.227 | 0.88 / 0.87 | 0.774 / 0.778 | 0.67 / 0.64 | 0.55 / 0.53 | -0.29 / -0.30 | -0.50 / -0.52 | 22.9 / 21.6 |
| Sc4_X5_d2 | -0.026 / -0.020 | 0.206 / 0.206 | 0.90 / 0.91 | 0.815 / 0.820 | 0.65 / 0.66 | 0.06 / 0.05 | -0.26 / -0.24 | -1.68 / -1.71 | 1.5 / 0.9 |
| Sc4_X7_d1 | -0.065 / -0.063 | 0.229 / 0.228 | 0.90 / 0.90 | 0.799 / 0.795 | 0.64 / 0.64 | 0.54 / 0.54 | -0.31 / -0.31 | -0.50 / -0.50 | 23.1 / 21.9 |
| Sc4_X7_d2 | -0.020 / -0.019 | 0.204 / 0.203 | 0.93 / 0.93 | 0.845 / 0.849 | 0.66 / 0.64 | 0.06 / 0.05 | -0.26 / -0.26 | -1.69 / -1.70 | 1.5 / 0.7 |
| Sc5_d1 | -0.029 / -0.025 | 0.212 / 0.210 | 0.93 / 0.90 | | | 0.12 / 0.10 | | -0.96 / -0.97 | 1.5 / 0.9 |
| Sc5_d2 | -0.027 / -0.019 | 0.203 / 0.205 | 0.93 / 0.91 | | | 0.01 / 0.01 | | -1.96 / -1.98 | 0.5 / 0.2 |

Paired v2 minus v1 differences with SEs (`per_replicate.csv`, 200 pairs). Real: the Sc3 bias, -0.013 to -0.014 (SE 0.002) at all three rho; the Sc4 delta 2 g_Rc, -0.028 (SE 0.011), and the Sc5 g, -0.012 to -0.014 (SE 0.003), small movements toward the truth; the Sc4 delta 1 bias, +0.007 (SE 0.003). Not real: the Sc1 squared error, +0.0016 (SE 0.0009), the Sc4 delta 2 map in R, +0.006 (SE 0.010), the Sc4 delta 2 g_R, +0.018 (SE 0.010), and every surface RMSE (within 0.010, SE 0.004 to 0.007). The v1 configuration is the Hg10 control: at R = 200 it is v1 itself; in the Phase 2b confirmation at R = 60 the two configurations differ by Sc3 bias 0.060 against 0.044, Sc4 delta 1 bias -0.053 against -0.040, Sc4 delta 2 map_R 0.67 against 0.69, all within one to two SEs, consistent with the R = 200 paired result. The tuning bought the Sc3 bias and nothing else; Hg10 is at least as good in Sc1 and at delta 1.

## Fractional ESS ladder (f50, f90; f25 not run)

The absolute ladder 50 / 75 / 100 and the fractional 0.5 / 0.9 of the ceiling give the same results wherever both are feasible: in Sc1 RMSE 0.182 / 0.179 / 0.175 against 0.179 / 0.172, with ESS_0 42.6 / 62.0 / 81.8 against 51.3 / 87.0, and f90 never caps where LRC-BART-100 caps in 24 percent; in Sc3, Sc4 and Sc5 every rung is within 0.005 in RMSE and 0.01 in the maps of LRC-BART-100. Sc2 is the case for the fractional ladder: the absolute rungs collapse to the capped 24.0 with realized ESS 13.2 to 13.7 and RMSE 0.183 to 0.186, while f50 and f90 give ESS_0 9.4 and 15.8, realized 5.9 and 9.5, RMSE 0.190 and 0.185, so only fractions separate the rungs there. In Sc4 at delta 1 the lower rungs have less bias (-0.055 f50, -0.056 at 50, against -0.064) and lower RMSE (0.219, 0.220 against 0.227), since less prior weight on the RWD leaves less moderate-shift leak.

## Why the boundary leak is not tunable

NOTES.md ("Phase 2b tuning") varied H_g from 1 to 30, the split prior, a slab up to 16 times the range rule, the leaf minimum, five sweeps per iteration, a warm start from a greedy tree and a 5000-iteration burn-in, and every configuration left g_R between -0.24 and -0.31 and g_Rc between -1.67 and -1.76. With f known (oracle), the g ensemble alone gives -0.13 to -0.19 in R and -1.86 to -1.93 in R^c, and by distance from the boundary the g error is -0.71 within 0.15 of X5 = 2, -0.37 at 0.15 to 0.3 and -0.17 beyond 0.5. The leak is boundary-local: the g cut is learned from about four RCT controls per 0.1 of X5, and f, fit to pooled data, smears the RWD step within 0.3 of the boundary, which g then compensates. Moderate shifts (Sc3, Sc4 delta 1) are a prior trade-off: the split-and-slab state pays the tree-rule prior (about 6.9 nats) plus the slab's Occam cost against a likelihood gain of order n delta^2 / (2 sigma^2), so the posterior mixes the detected state with states in which f absorbs part of the shift; fewer trees raise tau_1 and the slab cost, and a flat w prior lowers it, which is the only lever that moved Sc3. The manuscript should state that the compatibility map separates the regions (0.66 in R, 0.05 in R^c at delta 2; 0.85 / 0.84 at delta 0.5, where T2 predicts partial detection), that the surface and g criteria measure the f/g decomposition within 0.3 of a boundary with 100 RCT controls and are not met by any g configuration, that the surface comparison against BART-PP is reported both over the region and away from the boundary, and that the Sc3 bias of 0.05 and the Sc4 delta 1 loss of 0.024 in RMSE against BART-PP are the price of the spike.

## Verdict reasoning

The tuned configuration fixes the one failure that tuning could fix, the Sc3 bias, and the Sc1, Sc2, Sc3 and Sc5 results are those of a method that borrows when the sources agree and abstains when they do not, with C-BART, BART-CP and MAP failing where LRC-BART holds. The Sc4 criteria as written are not met and, on the evidence in NOTES.md, cannot be met by tuning g, so they will not be met at R = 500 either. The honest position for the paper is that local robustness is demonstrated on the map and the ATE, not on the control surface in R, and that against BART-PP the method ties at delta 2 and loses at delta 1; the manuscript must say this, and the full run must carry Hg10 and f25 and the paper's power definition so the sensitivity and the Sc1 power criterion are reportable. Proceed to the full Gaussian and survival study with H_g = 5, (0.5, 3), range-rule slab, Beta(1, 1) w, and the fractional ladder as the primary sensitivity. Criteria are left as written and reported as failed.
