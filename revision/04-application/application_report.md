# Phase 4 application: LRC-BART on the EloKRd + UCMM data

Date: 2026-09-09. Code in `code/`, summaries in `results/`, table in `results_table.csv`, figures in `figures/`. Two cores, about 15 minutes of compute.

## Data check

The n = 283 dataset reproduces every count and survival summary of Table baseline (30 and 253 patients, 13 and 10 EloKRd events, 5-year RMST 3.71 vs 3.65 for PFS and 4.12 vs 4.57 for OS, log-rank p = 0.506 and 0.048, KM ratio 1.015); see `data_check.md`. Two coding problems sit behind the matching counts. Ethnicity is labelled "Non-Hispanic" in EloKRd and "Not Hispanic or Latino" in UCMM, so Yunxuan's `prep_covariates()` makes two non-Hispanic dummies, one per source: in every model-based row the 29 non-Hispanic EloKRd patients are predicted from the 13-patient UCMM Hispanic stratum, and the Stan AFT rows add a prior-only coefficient on top. Race keeps Asian as a separate reference level although the table pools it into Other/Unknown. We run everything under the original coding (to reproduce the paper) and under a harmonized seven-column coding (age, male, Black, Other, Hispanic, high-risk cytogenetics, ASCT). The n78 file is the secondary KRd/Rd-only cohort and is not used.

## Reference rows

Under the original coding Yunxuan's scripts (KM, complete-pooling AFT, hierarchical AFT in Stan) and lrcbart's plain AFT-BART in place of the missing `cabart.cpp` reproduce Table rmst-results to rounding: PFS 1.015, 1.03, 1.24 and 1.06 against the published 1.02, 1.03, 1.21 and 1.06, with matching intervals; OS 0.90, 0.90, 1.25 (0.77, 2.34) and 0.91 against 0.90, 0.90, 1.18 (0.77, 1.89) and 0.90. The published hierarchical row is the s^2 = 0.5 prior, and the published point estimates are posterior medians. The MAP-AFT-BART rows cannot be rerun because `cmbart.cpp` is not in the repo.

Under the harmonized coding the complete-pooling AFT falls to 0.87 (0.69, 1.06) for PFS and 0.80 (0.65, 0.94) for OS, the only interval in the table that excludes one; standard BART falls to 0.97 (0.79, 1.16) and 0.905 (0.78, 1.02); the hierarchical AFT stays near 1.2. Part of the published upward drift of the parametric rows is the ethnicity coding.

## LRC-BART results

Single-arm mode: f is AFT-BART on UCMM (H_f = 10 as in the paper; H_f = 50 also run), g a prior draw with H_g = 5, split prior (0.5, 3), nu_0 = 3; the treated arm is AFT-BART with 50 trees; w = 1 reported, w = 0.9 sensitivity; 4 chains of 2000 warm-up and 2000 draws. Calibration on UCMM with sigma_1 = sigma_2 gives a ceiling sigma_1^2 / V_mu^f of 63.7 for PFS and 36.7 for OS, far below 253, because the 30 profiles (30 percent ASCT, 50 percent high-risk) lie where UCMM is thin. The absolute ladder is feasible: s_0^2 = 0.0056, 0.010, 0.018 for PFS (ESS_0 = 28.1, 20.6, 14.2) and 0.0050, 0.010, 0.022 for OS (24.4, 19.0, 12.7); the fractional ladder (0.9 / 0.5 / 0.25 of the ceiling, all rows in `results_table.csv`) spans N_target 57 to 16 (PFS) and 33 to 9 (OS). Under the original coding the ceiling collapses to 12.6 and 14.6 and N_target = 30 is capped, so the calibration flags the coding problem on its own.

Standardized 5-year RMST ratio, posterior mean (95 percent CrI), harmonized coding, H_f = 10, w = 1:

| Target | PFS | OS |
|---|---|---|
| N_target = 30 | 0.974 (0.774, 1.227) | 0.903 (0.765, 1.043) |
| N_target = 22.5 | 0.980 (0.764, 1.294) | 0.905 (0.763, 1.064) |
| N_target = 15 | 0.989 (0.751, 1.401) | 0.911 (0.757, 1.108) |
| 0.9 x ceiling | 0.967 (0.791, 1.153) | 0.902 (0.766, 1.036) |

P(ratio > 1) is 0.37 (PFS) and 0.08 (OS) at N_target = 30. The median ratio is 1.10 (0.49, 2.61) for PFS and 0.55 (0.22, 1.40) for OS, wide because neither EloKRd median is reached. Sensitivities: w = 0.9 widens PFS to 0.997 (0.741, 1.477) and OS to 0.909 (0.764, 1.098), the slab adding prior noise as Section F predicts; H_f = 50 changes means by at most 0.014; the original coding gives 1.04 (0.79, 1.43) and 0.91 (0.76, 1.07). Realized ESS equals ESS_0 within one unit in every row. Convergence: split-Rhat 1.000 for every estimand, bulk ESS 6800 to 7900 and tail ESS 5300 to 7300 of 8000 draws for the RMST ratio, bulk ESS at least 2400 for the median ratio.

Figures. `fig_ess_map`: ESS(x) at the 30 profiles is 5.3 to 16.9 (median 9.5) for PFS and 5.8 to 14.9 for OS, highest for White patients with ASCT and standard-risk cytogenetics, lowest for the one Hispanic patient, the Other-race patients and the high-risk, no-ASCT profiles. `fig_rmst_posterior`: the three posteriors of the RMST ratio. `fig_g_profiles`: g(x) at the 30 profiles is a prior draw in single-arm mode, mean zero everywhere, SD 0.26 on the log-time scale at N_target = 30 (0.10 at 0.9 x ceiling, 0.43 at N_target = 15).

## Comparison with the paper's MAP-AFT-BART rows

The long upper tails do not persist. The paper's 97.5 percent limits are 16.4, 21.8, 32.8 for PFS and 4.77, 5.55, 7.22 for OS; ours are 1.23, 1.29, 1.40 and 1.04, 1.06, 1.11. Over 8000 draws the PFS ratio at N_target = 30 has a 99.9 percentile of 1.61, a maximum of 3.09 and P(ratio > 2) = 0.0003 (0.002 at N_target = 15); no OS draw exceeds 2.04. The paper's tails came from the per-leaf discrepancy prior letting the control RMST approach zero on rare draws; with the discrepancy carried by a five-tree g ensemble of variance H_g tau_0^2 (0.07 on the log scale at N_target = 30) that route is closed. Point estimates move from 1.08 to 1.10 down to 0.97 to 0.99 for PFS and from 0.97 to 1.00 down to 0.90 to 0.92 for OS; about 0.06 of the PFS change is the ethnicity coding (the original-coding LRC-BART run gives 1.04), the rest the model. Interval width still grows as the target falls, so the calibration trade-off is intact, on a scale of tenths rather than tens.

## What Section 5 should say

State the harmonization fix and rerun the reference rows under it, noting that the complete-pooling AFT interval for OS then excludes one while every nonparametric row includes it. Report the ESS ceiling (64 for PFS, 37 for OS) first: it is the amount of UCMM information that reaches the 30 profiles, and it explains why a target of 30 is worth about half the usable registry. Replace the MAP-AFT-BART rows by the LRC-BART rows at N_target = 30, 22.5 and 15, with the fractional ladder in the supplement, report posterior means with medians alongside, and keep the conclusion of no evidence of a marginal EloKRd benefit, now with upper limits below 1.5. Add the local ESS map as the covariate-resolved figure, state that g is a prior draw here so the borrowing map has no single-arm counterpart, and drop the sentence attributing the long upper tails to an identification limit: they were a property of the earlier prior, not of the data.
