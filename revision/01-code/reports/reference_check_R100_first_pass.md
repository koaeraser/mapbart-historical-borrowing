# Reference check: DGP and reference comparators against the paper's Table 1

Run: 2026-09-09 11:03. Gaussian outcome, scenarios Sc1, Sc2, Sc3 (rho = 0), R = 100 replicates, seed 2026, 18 workers, 1.3 minutes wall time. Paper values are from paper/tables_figures.tex (Table gaussian_results, R = 500). Tolerances requested: bias 0.03, RMSE 0.03, coverage 0.03. RMSE here is the Monte Carlo RMSE of the posterior mean (see R/metrics.R for why this, and not the draw-based formula of the paper's metric table, is what the published numbers correspond to).

Verdict: FAIL. 14 of 54 checked cells (bias, RMSE, coverage for 18 scenario x method rows) fall outside the tolerance.

## ATE metrics (paper / ours)

| Scenario | Method | Bias paper | Bias ours | SE | SD paper | SD ours | RMSE paper | RMSE ours | Cov paper | Cov ours | Pow paper | Pow ours | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | -0.03 | 0.01 | 0.036 | 0.34 | 0.34 | 0.35 | 0.36 | 0.95 | 0.93 | 0.41 | 0.81 | bias |
| Sc1 | LM-CP | -0.02 | 0.03 | 0.025 | 0.24 | 0.24 | 0.26 | 0.25 | 0.93 | 0.94 | 0.62 | 0.99 | bias |
| Sc1 | LM-PP | -0.02 | -0.00 | 0.036 | 0.34 | 0.35 | 0.34 | 0.35 | 0.96 | 0.93 | 0.40 | 0.82 | cov |
| Sc1 | BART-NP | -0.01 | 0.00 | 0.022 | 0.20 | 0.20 | 0.22 | 0.21 | 0.92 | 0.92 | 0.77 | 1.00 | ok |
| Sc1 | BART-CP | 0.00 | 0.00 | 0.014 | 0.14 | 0.14 | 0.15 | 0.14 | 0.92 | 0.98 | 0.95 | 1.00 | cov |
| Sc1 | BART-PP | 0.00 | 0.00 | 0.019 | 0.21 | 0.21 | 0.20 | 0.19 | 0.96 | 0.96 | 0.78 | 1.00 | ok |
| Sc2 | LM-NP | -0.04 | 0.01 | 0.034 | 0.32 | 0.32 | 0.33 | 0.34 | 0.94 | 0.94 | 0.44 | 0.86 | bias |
| Sc2 | LM-CP | 0.37 | 0.40 | 0.031 | 0.28 | 0.29 | 0.47 | 0.50 | 0.73 | 0.72 | 0.91 | 1.00 | rmse |
| Sc2 | LM-PP | 0.02 | 0.03 | 0.036 | 0.33 | 0.34 | 0.33 | 0.36 | 0.96 | 0.92 | 0.48 | 0.83 | cov |
| Sc2 | BART-NP | -0.02 | 0.01 | 0.023 | 0.20 | 0.20 | 0.21 | 0.23 | 0.93 | 0.94 | 0.77 | 1.00 | ok |
| Sc2 | BART-CP | 0.05 | 0.06 | 0.019 | 0.18 | 0.18 | 0.18 | 0.20 | 0.95 | 0.88 | 0.93 | 1.00 | cov |
| Sc2 | BART-PP | 0.00 | 0.02 | 0.021 | 0.21 | 0.21 | 0.19 | 0.21 | 0.96 | 0.97 | 0.79 | 1.00 | ok |
| Sc3_rho0 | LM-NP | -0.03 | 0.04 | 0.034 | 0.36 | 0.37 | 0.37 | 0.34 | 0.94 | 0.99 | 0.39 | 0.83 | bias,rmse,cov |
| Sc3_rho0 | LM-CP | 1.01 | 1.04 | 0.024 | 0.26 | 0.26 | 1.04 | 1.07 | 0.02 | 0.02 | 1.00 | 1.00 | ok |
| Sc3_rho0 | LM-PP | 0.06 | 0.07 | 0.033 | 0.35 | 0.36 | 0.36 | 0.34 | 0.94 | 0.98 | 0.49 | 0.87 | cov |
| Sc3_rho0 | BART-NP | -0.01 | 0.04 | 0.029 | 0.29 | 0.29 | 0.31 | 0.29 | 0.93 | 0.93 | 0.52 | 0.93 | bias |
| Sc3_rho0 | BART-CP | 0.97 | 0.99 | 0.020 | 0.20 | 0.20 | 1.00 | 1.01 | 0.01 | 0.00 | 1.00 | 1.00 | ok |
| Sc3_rho0 | BART-PP | 0.02 | 0.06 | 0.027 | 0.28 | 0.28 | 0.30 | 0.28 | 0.93 | 0.95 | 0.58 | 0.97 | bias |

## Arm-specific metrics (paper / ours)

| Scenario | Method | Trt bias paper | ours | Trt RMSE paper | ours | Ctrl bias paper | ours | SE | Ctrl RMSE paper | ours |
|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | -0.01 | -0.01 | 0.24 | 0.16 | 0.02 | -0.01 | 0.025 | 0.36 | 0.25 |
| Sc1 | LM-CP | -0.01 | -0.01 | 0.24 | 0.16 | 0.01 | -0.03 | 0.018 | 0.22 | 0.18 |
| Sc1 | LM-PP | -0.01 | -0.01 | 0.24 | 0.16 | 0.00 | -0.01 | 0.025 | 0.36 | 0.25 |
| Sc1 | BART-NP | 0.00 | -0.00 | 0.15 | 0.12 | 0.01 | -0.01 | 0.017 | 0.24 | 0.17 |
| Sc1 | BART-CP | 0.00 | 0.00 | 0.15 | 0.12 | 0.00 | -0.00 | 0.009 | 0.13 | 0.09 |
| Sc1 | BART-PP | 0.00 | -0.00 | 0.15 | 0.12 | 0.00 | -0.01 | 0.015 | 0.23 | 0.15 |
| Sc2 | LM-NP | -0.01 | 0.01 | 0.22 | 0.14 | 0.03 | -0.00 | 0.026 | 0.34 | 0.26 |
| Sc2 | LM-CP | -0.01 | 0.01 | 0.22 | 0.14 | -0.38 | -0.38 | 0.024 | 0.46 | 0.45 |
| Sc2 | LM-PP | -0.01 | 0.01 | 0.22 | 0.14 | -0.03 | -0.02 | 0.027 | 0.36 | 0.27 |
| Sc2 | BART-NP | 0.00 | 0.01 | 0.15 | 0.11 | 0.02 | -0.00 | 0.020 | 0.23 | 0.20 |
| Sc2 | BART-CP | 0.00 | 0.01 | 0.15 | 0.11 | -0.05 | -0.06 | 0.016 | 0.19 | 0.17 |
| Sc2 | BART-PP | 0.00 | 0.01 | 0.15 | 0.11 | 0.00 | -0.02 | 0.017 | 0.23 | 0.17 |
| Sc3_rho0 | LM-NP | 0.00 | 0.01 | 0.25 | 0.15 | 0.03 | -0.02 | 0.024 | 0.38 | 0.24 |
| Sc3_rho0 | LM-CP | 0.00 | 0.01 | 0.25 | 0.15 | -1.01 | -1.02 | 0.018 | 1.02 | 1.04 |
| Sc3_rho0 | LM-PP | 0.00 | 0.01 | 0.25 | 0.15 | -0.06 | -0.05 | 0.024 | 0.37 | 0.24 |
| Sc3_rho0 | BART-NP | 0.02 | 0.03 | 0.20 | 0.14 | 0.03 | -0.01 | 0.022 | 0.33 | 0.22 |
| Sc3_rho0 | BART-CP | 0.02 | 0.02 | 0.20 | 0.14 | -0.96 | -0.97 | 0.014 | 0.96 | 0.98 |
| Sc3_rho0 | BART-PP | 0.02 | 0.03 | 0.20 | 0.14 | 0.00 | -0.03 | 0.021 | 0.30 | 0.21 |

## Largest discrepancies

| Scenario | Method | d bias | d SD | d RMSE | d cov | d pow | d ctrl bias | d ctrl RMSE | d trt RMSE |
|---|---|---|---|---|---|---|---|---|---|
| Sc2 | BART-CP | 0.01 | -0.00 | 0.02 | -0.07 | 0.07 | -0.01 | -0.02 | -0.04 |
| Sc3_rho0 | LM-NP | 0.07 | 0.01 | -0.03 | 0.05 | 0.44 | -0.05 | -0.14 | -0.10 |
| Sc1 | BART-CP | 0.00 | 0.00 | -0.01 | 0.06 | 0.05 | -0.00 | -0.04 | -0.03 |
| Sc2 | LM-NP | 0.05 | -0.00 | 0.01 | 0.00 | 0.42 | -0.03 | -0.08 | -0.08 |
| Sc1 | LM-CP | 0.05 | 0.00 | -0.01 | 0.01 | 0.37 | -0.04 | -0.04 | -0.08 |

## Notes

SE is the Monte Carlo standard error of our bias estimate (R replicates); the paper's R = 500 values carry roughly half that SE. A cell flagged OUT whose discrepancy is within about two combined standard errors is consistent with Monte Carlo variation. Systematic discrepancies (the same sign across all scenarios, or larger than three SEs) point to a difference in the mechanism or the estimator and are discussed in the diagnosis section below, written by hand after inspecting this table.

## Diagnosis

(to be completed after inspection)
