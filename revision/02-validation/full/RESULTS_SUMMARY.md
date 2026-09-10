# Phase 2 full study: results summary

Date 2026-09-09. R = 500, seed 2026 (replicates 1 to 200 paired with the dev run), 16 workers, configuration H_g = 5, split prior (0.5, 3), range-rule slab, w ~ Beta(1, 1). No fit errors in any block. Wall time: Gaussian 178.7 min, survival 118.3 min, single-arm 36.8 min, total 333.8 min (5.6 h). Files: `<block>/summary_table.csv`, `<block>/wall_time.txt`, the full tables and change lists in `analysis_all.md` (from `analyze_full.R`). `pow` is the paper's stated definition (95% interval excludes the null), `pow90` the 90% interval; in the Gaussian study both equal 1.00 for every BART method in every scenario and separate only the linear methods (LM-PP Sc3 0.85 against 0.91), so the published power column (0.77 for BART-PP) is not reproducible under either definition at delta = 1.

## Gaussian (12 scenarios, 20 methods)

Bias, RMSE and coverage of the ATE; best linear method in the last row of each block.

| scenario | LRC-BART-100 | BART-PP | BART-NP | BART-CP | C-BART | best linear |
|---|---|---|---|---|---|---|
| Sc1 | 0.009 / 0.172 / 0.94 | 0.006 / 0.205 / 0.94 | 0.004 / 0.220 / 0.91 | 0.008 / 0.142 / 0.96 | 0.007 / 0.150 / 0.95 | LM-CP 0.009 / 0.236 / 0.96 |
| Sc2 | 0.023 / 0.192 / 0.94 | 0.007 / 0.204 / 0.96 | -0.008 / 0.219 / 0.93 | 0.058 / 0.195 / 0.91 | 0.042 / 0.186 / 0.91 | PSCL 0.024 / 0.292 / 0.97 |
| Sc3 rho 0 | 0.041 / 0.280 / 0.94 | 0.029 / 0.273 / 0.94 | -0.002 / 0.289 / 0.94 | 0.984 / 1.003 / 0.00 | 0.226 / 0.380 / 0.87 | LM-PP 0.043 / 0.346 / 0.95 |
| Sc3 rho +-0.5 | 0.040 to 0.042 / 0.275 to 0.278 / 0.93 | 0.027 to 0.029 / 0.267 to 0.270 / 0.93 | -0.004 to 0.001 / 0.280 to 0.282 / 0.94 | 0.846 / 0.871 / 0.01 | 0.28 / 0.43 / 0.81 | LM-PP 0.041 / 0.342 / 0.95 |
| Sc4 X5 d0.5 | -0.043 / 0.200 / 0.92 | 0.000 / 0.203 / 0.95 | 0.004 / 0.220 / 0.91 | -0.164 / 0.217 / 0.78 | -0.109 / 0.190 / 0.90 | PSCL -0.120 / 0.287 / 0.95 |
| Sc4 X5 d1 | -0.041 / 0.223 / 0.89 | -0.002 / 0.207 / 0.95 | 0.004 / 0.220 / 0.91 | -0.338 / 0.368 / 0.35 | -0.207 / 0.273 / 0.73 | LM-PP -0.012 / 0.339 / 0.95 |
| Sc4 X5 d2 | 0.000 / 0.202 / 0.92 | -0.006 / 0.205 / 0.96 | 0.004 / 0.220 / 0.91 | -0.685 / 0.701 / 0.00 | -0.227 / 0.391 / 0.66 | LM-NP -0.002 / 0.342 / 0.94 |
| Sc4 X7 d1, d2 | -0.038 / 0.224 / 0.90; 0.001 / 0.203 / 0.93 | -0.001 / 0.206 / 0.96; -0.002 / 0.206 / 0.95 | as Sc1 | -0.340 / 0.369; -0.685 / 0.701 | -0.206 / 0.276; -0.218 / 0.383 | LM-PP -0.012 / 0.337; -0.023 / 0.339 |
| Sc5 d1, d2 | -0.005 / 0.210 / 0.91; -0.001 / 0.203 / 0.92 | -0.004 / 0.205 / 0.95; -0.012 / 0.205 / 0.95 | as Sc1 | -0.687 / 0.702; -1.390 / 1.398 | -0.244 / 0.373; -0.316 / 0.650 | LM-PP -0.022 / 0.338; -0.042 / 0.339 |

Criteria: 18 of 25 pass, the same count as the dev run, with two flips. Item 1 (Sc3 rho -0.5 bias) moves from 0.050 to 0.042 and passes; item 18 (mean g in R within 0.25 of 0 at delta 2) moves from -0.244 to -0.255 and fails. The seven failures are the Sc1 power criterion (both at 1.00), the five Sc4 delta 2 local-borrowing criteria (surface RMSE in R 0.824 against 0.9 x 0.785; all-spike map in R 0.269; P(|g| < 0.5) in R 0.643; g in R -0.255, in R^c -1.713) and the Sc4 delta 1 RMSE (0.223 against BART-PP 0.207, coverage 0.89). At R = 500 the Sc3 bias of LRC-BART-100 is 0.040 to 0.042 against BART-PP's 0.027 to 0.029 and BART-NP's near zero, with RMSE within 0.008 of BART-PP.

Sc4 map and surface. The compatibility map separates the regions at delta 2 (0.64 in R, 0.05 in R^c, both X5 and X7) and partially at delta 1 (0.62 / 0.52) and delta 0.5 (0.84 / 0.83); Sc5 abstains (0.10 at delta 1, 0.01 at delta 2). The control-surface RMSE at delta 2 is 0.824 in R and 0.857 in R^c against BART-PP's 0.785 and 0.822 (X5), 0.850 / 0.848 against 0.805 / 0.806 (X7); at delta 0.5 and 1 LRC-BART is within 0.02 of BART-PP in R (0.751 vs 0.758, 0.783 vs 0.784) and behind it in R^c by 0.001 and 0.023. C-BART's surface in R^c is 1.087 at delta 2. Hg10 (the dev-run configuration) is within 0.005 of Hg5 on every Sc4 quantity and worse on the Sc3 bias (0.053 to 0.055).

Fractional ladder. In Sc1 the rungs order as expected (RMSE 0.184 / 0.178 / 0.171 for f25 / f50 / f90, 0.181 / 0.176 / 0.172 for 50 / 75 / 100; prior ESS 28 / 54 / 90 against 43 / 62 / 82, ceiling 127) and LRC-BART-100 is capped in 20% of replicates while f90 never is. In Sc2 the absolute rungs collapse to the ceiling (ESS_0 24.5, capped 100%, bias 0.022 to 0.023) while f25 / f50 / f90 give ESS_0 5.1 / 9.7 / 16.5 with bias 0.012 / 0.016 / 0.019 and RMSE 0.200 / 0.197 / 0.193. Under the moderate Sc4 shift the lower rungs are less biased (f25 bias -0.032, RMSE 0.208 at delta 1 against -0.041, 0.223 at 100). In Sc3 and Sc5 every rung is within 0.005 in bias and RMSE of LRC-BART-100.

## Survival (8 scenarios, 15 methods)

Median ratio (bias / RMSE / coverage); RMST ratio at tau = 3 in parentheses for LRC-BART-100 and BART-PP.

| scenario | LRC-BART-100 | BART-PP | BART-NP | BART-CP | C-BART | best AFT |
|---|---|---|---|---|---|---|
| Sc1 | 0.045 / 0.312 / 0.95 (0.003 / 0.097) | 0.032 / 0.345 / 0.95 (-0.023 / 0.106) | 0.260 / 0.511 / 0.92 | 0.064 / 0.232 / 0.96 | 0.034 / 0.248 / 0.97 | AFT-CP 0.165 / 0.428 / 0.95 |
| Sc2 | 0.061 / 0.322 / 0.94 (0.018 / 0.100) | 0.036 / 0.335 / 0.95 (-0.002 / 0.097) | 0.251 / 0.495 / 0.91 | 0.194 / 0.383 / 0.93 | 0.138 / 0.331 / 0.93 | AFT-NP 0.175 / 0.672 / 0.96 |
| Sc3 (three rho) | 0.147 to 0.155 / 0.539 to 0.558 / 0.93 to 0.95 (0.030 to 0.032 / 0.144 to 0.146) | 0.041 to 0.047 / 0.441 to 0.451 / 0.91 to 0.92 (0.000 to -0.005 / 0.134 to 0.137) | 0.186 to 0.195 / 0.555 to 0.574 / 0.90 | 1.95 to 2.31 / 2.06 to 2.41 / 0.00 | 0.58 to 0.63 / 0.91 to 0.95 / 0.84 to 0.88 | AFT-NP 0.22 to 0.23 / 0.81 to 0.83 / 0.94 |
| Sc4 X5 d0.5 | 0.001 / 0.317 / 0.95 | 0.011 / 0.339 / 0.96 | 0.260 / 0.511 | -0.168 / 0.252 / 0.86 | -0.094 / 0.247 / 0.94 | AFT-CP -0.123 / 0.347 / 0.92 |
| Sc4 X5 d1 | -0.011 / 0.340 / 0.93 (-0.037 / 0.108) | -0.012 / 0.339 / 0.97 (-0.059 / 0.113) | 0.260 / 0.511 | -0.367 / 0.401 / 0.50 | -0.202 / 0.305 / 0.84 | AFT-PP (RMST) -0.039 / 0.139 |
| Sc5 d1 | 0.032 / 0.353 / 0.93 (0.001 / 0.112) | -0.029 / 0.333 / 0.95 (-0.041 / 0.110) | 0.260 / 0.511 | -0.703 / 0.712 / 0.00 | -0.233 / 0.405 / 0.77 | AFT-PP 0.193 / 0.695 / 0.94 |

The median-ratio Sc3 bias of LRC-BART (0.15) is three times BART-PP's, with coverage at nominal; on the RMST scale it is 0.03 against 0.00, RMSE 0.145 against 0.135. In Sc4 (log shift 1, residual scale 1.4) the map separates weakly (0.64 in R, 0.58 in R^c) and the surface in R is 0.792 against BART-PP's 0.803, in R^c 0.868 against 0.888. Survival power is low for every method (LRC-BART 0.34 in Sc1, BART-PP 0.25; pow90 about 0.1 higher). One comparator differs from the paper's Table 2 and should be checked before writing: BART-NP's median-ratio bias is 0.26 (paper 0.01) with posterior SD 0.42 (paper 0.16); the survival comparators were never reference-checked. The fractional ladder behaves as in the Gaussian study (Sc2 rungs collapse at the ceiling of 40; f25 bias 0.042 against 0.061).

## Single-arm (paper's Tables 5 to 7)

Bias / RMSE / coverage; the paper's value in brackets.

| design | scenario | BART-CP | LRC-BART-full (paper: near-complete pooling) | LRC-BART-f90-w1 | LM-CP or AFT-CP |
|---|---|---|---|---|---|
| Gaussian n 200 | Sc1 | 0.010 / 0.162 / 0.94 [-0.01 / 0.16 / 0.93] | 0.007 / 0.162 / 0.94 [-0.01 / 0.16 / 0.97] | 0.007 / 0.163 / 0.98 | 0.007 / 0.254 / 0.95 [-0.02 / 0.25 / 0.95] |
| Gaussian n 200 | Sc2 | 0.408 / 0.506 / 0.75 [0.05 / 0.19 / 0.92] | 0.420 / 0.515 / 0.75 [0.05 / 0.19 / 0.95] | 0.417 / 0.511 / 0.89 | 1.682 / 1.734 / 0.01 [0.28 / 0.41 / 0.79] |
| Survival n 30, median | Sc1 | -0.013 / 0.482 / 0.97 [0.04 / 0.47 / 1.00] | 0.182 / 0.686 / 0.98 [0.04 / 0.47 / 1.00] | 0.222 / 0.711 / 0.99 | 1.540 / 4.444 / 0.92 [0.11 / 0.62 / 0.97] |
| Survival n 30, median | Sc2 | 0.814 / 1.231 / 0.92 [0.18 / 0.58 / 0.99] | 1.324 / 1.864 / 0.95 [0.19 / 0.59 / 0.99] | 1.511 / 2.055 / 0.97 | 20.9 / 51.6 / 0.19 [0.80 / 1.37 / 0.87] |
| Survival n 200, median | Sc1 | 0.108 / 0.286 / 0.95 [-0.00 / 0.20 / 0.96] | 0.076 / 0.262 / 0.97 [-0.00 / 0.21 / 0.99] | 0.087 / 0.268 / 0.99 | 0.168 / 0.480 / 0.95 [0.07 / 0.38 / 0.96] |
| Survival n 200, median | Sc2 | 1.136 / 1.359 / 0.60 [0.08 / 0.25 / 0.95] | 1.195 / 1.409 / 0.64 [0.08 / 0.25 / 0.98] | 1.336 / 1.562 / 0.82 | 9.478 / 10.85 / 0.01 [0.71 / 0.95 / 0.76] |
| Survival n 200, RMST | Sc1 | -0.001 / 0.082 / 0.97 [-0.01 / 0.07 / 0.96] | -0.024 / 0.076 / 0.95 [-0.01 / 0.07 / 0.98] | -0.022 / 0.076 / 0.99 | -0.023 / 0.098 / 0.94 [-0.01 / 0.10 / 0.96] |
| Survival n 200, RMST | Sc2 | 0.305 / 0.361 / 0.61 [0.04 / 0.09 / 0.93] | 0.284 / 0.336 / 0.70 [0.04 / 0.09 / 0.98] | 0.301 / 0.353 / 0.85 | 1.164 / 1.244 / 0.01 [0.17 / 0.22 / 0.71] |

Sc1 reproduces the paper for BART-CP, LM-CP and the near-complete-pooling model at n 200. Sc2 does not: with the paper's stated mechanism (Sc2 source assignment, trial treatment-only, RWD the sole control) every method is far more biased than the published table (LM-CP 1.68 against 0.28, BART-CP 0.41 against 0.05), while the two-arm Sc2 reproduced Table 1 (LM-CP 0.39 against 0.37). The published single-arm Sc2 rows were therefore not produced by the stated design; the manuscript must state the mechanism that produced them or replace them. LRC-BART-100-w1 equals LRC-BART-full in every design because N_target = 100 sits at or above the single-arm ESS ceiling (Sc2 ceiling 33 Gaussian, 27 to 34 survival) or, in Sc1 at n 200, within 0.005 of it. The w = 0.9 sensitivity reproduces the paper's "N_target = 100" pathology (posterior SD 1.38 Gaussian, 10 to 34 on the median ratio, power 0, coverage 1.00) with the point estimate unchanged, which supports w = 1 as the reported single-arm model. HierLM / HierAFT with the paper's stated hyperprior (tau_j^2 ~ IG(1.5, 0.75)) and no RCT controls draw the control surface from the hierarchy: SD 8.2 for the Gaussian ATE and median-ratio means of order 10^50; the published SD of 0.58 and 1.4 to 2.7 cannot come from that hyperprior, so these rows are reported as unusable, not tuned to match.

## Changes relative to the dev run (R = 200), Gaussian

Criterion flips: item 1 FAIL to PASS (Sc3 rho -0.5 bias 0.050 to 0.042), item 18 PASS to FAIL (g in R -0.244 to -0.255). Threshold crossings (values in `analysis_all.md`): |bias| falls below 0.05 for LRC-BART-100 / 75 / 50 / f50 in Sc3 and for every LRC rung in Sc4 X5 delta 0.5 and 1 and X7 delta 1 (dev -0.051 to -0.067, full -0.031 to -0.043), for HierLM in three scenarios and LM-PP in two; MAP-50 falls below 0.10 in Sc4 delta 0.5. Coverage rises above 0.92 for the lower LRC rungs in Sc4 delta 0.5 and 2 and Sc5, falls below it for C-BART in Sc2 (0.906) and w0 in Sc4 X7 delta 1; rises above 0.90 for LRC-BART-100 / f90 in Sc4 delta 0.5 (0.880 to 0.918 / 0.926) and BART-NP everywhere (0.912), falls below it for MAP-75 in Sc3 (0.884). RMSE order against BART-PP reverses in LRC-BART's favour for the 100 / 75 / f90 rungs and C-BART in Sc4 delta 0.5 and for four rungs at delta 2 (0.202 against 0.205). Sign changes are confined to biases below 0.03 in absolute value (Sc1 for every method, BART-NP and w0 in Sc4 and Sc5) and to control biases below 0.02. The largest headline movements are C-BART in Sc5 delta 2 (-0.350 to -0.316) and the LRC-BART moderate-shift biases in Sc4 (about +0.023).
