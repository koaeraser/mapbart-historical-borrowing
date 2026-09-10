# Reference check: DGP and reference comparators against the paper's Table 1

Run: 2026-09-09 11:10. Gaussian outcome, scenarios Sc1, Sc2, Sc3 (rho = 0), R = 500 replicates, seed 2026, 18 workers, 5.2 minutes wall time. Paper values are from paper/tables_figures.tex (Table gaussian_results, R = 500). Tolerances requested: bias 0.03, RMSE 0.03, coverage 0.03. RMSE here is the Monte Carlo RMSE of the posterior mean (see R/metrics.R for why this, and not the draw-based formula of the paper's metric table, is what the published numbers correspond to).

Verdict: PASS with exceptions. 3 of 54 checked cells (bias, RMSE, coverage for 18 scenario x method rows) fall outside the tolerance.

## ATE metrics (paper / ours)

| Scenario | Method | Bias paper | Bias ours | SE | SD paper | SD ours | RMSE paper | RMSE ours | Cov paper | Cov ours | Pow paper | Pow ours | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | -0.03 | -0.00 | 0.015 | 0.34 | 0.34 | 0.35 | 0.34 | 0.95 | 0.94 | 0.41 | 0.83 | ok |
| Sc1 | LM-CP | -0.02 | 0.01 | 0.011 | 0.24 | 0.24 | 0.26 | 0.24 | 0.93 | 0.96 | 0.62 | 0.99 | ok |
| Sc1 | LM-PP | -0.02 | -0.00 | 0.015 | 0.34 | 0.35 | 0.34 | 0.34 | 0.96 | 0.94 | 0.40 | 0.85 | ok |
| Sc1 | BART-NP | -0.01 | 0.00 | 0.010 | 0.20 | 0.20 | 0.22 | 0.22 | 0.92 | 0.91 | 0.77 | 1.00 | ok |
| Sc1 | BART-CP | 0.00 | 0.01 | 0.006 | 0.14 | 0.14 | 0.15 | 0.14 | 0.92 | 0.96 | 0.95 | 1.00 | cov |
| Sc1 | BART-PP | 0.00 | 0.01 | 0.009 | 0.21 | 0.21 | 0.20 | 0.20 | 0.96 | 0.94 | 0.78 | 1.00 | ok |
| Sc2 | LM-NP | -0.04 | -0.02 | 0.014 | 0.32 | 0.32 | 0.33 | 0.32 | 0.94 | 0.95 | 0.44 | 0.87 | ok |
| Sc2 | LM-CP | 0.37 | 0.38 | 0.014 | 0.28 | 0.29 | 0.47 | 0.49 | 0.73 | 0.73 | 0.91 | 0.99 | ok |
| Sc2 | LM-PP | 0.02 | 0.01 | 0.015 | 0.33 | 0.34 | 0.33 | 0.33 | 0.96 | 0.95 | 0.48 | 0.84 | ok |
| Sc2 | BART-NP | -0.02 | -0.01 | 0.010 | 0.20 | 0.20 | 0.21 | 0.22 | 0.93 | 0.93 | 0.77 | 1.00 | ok |
| Sc2 | BART-CP | 0.05 | 0.06 | 0.008 | 0.18 | 0.18 | 0.18 | 0.19 | 0.95 | 0.91 | 0.93 | 1.00 | cov |
| Sc2 | BART-PP | 0.00 | 0.01 | 0.009 | 0.21 | 0.21 | 0.19 | 0.20 | 0.96 | 0.96 | 0.79 | 1.00 | ok |
| Sc3_rho0 | LM-NP | -0.03 | 0.01 | 0.016 | 0.36 | 0.37 | 0.37 | 0.36 | 0.94 | 0.96 | 0.39 | 0.80 | bias |
| Sc3_rho0 | LM-CP | 1.01 | 1.04 | 0.011 | 0.26 | 0.26 | 1.04 | 1.07 | 0.02 | 0.03 | 1.00 | 1.00 | ok |
| Sc3_rho0 | LM-PP | 0.06 | 0.04 | 0.015 | 0.35 | 0.36 | 0.36 | 0.35 | 0.94 | 0.95 | 0.49 | 0.85 | ok |
| Sc3_rho0 | BART-NP | -0.01 | -0.00 | 0.013 | 0.29 | 0.29 | 0.31 | 0.29 | 0.93 | 0.94 | 0.52 | 0.92 | ok |
| Sc3_rho0 | BART-CP | 0.97 | 0.98 | 0.009 | 0.20 | 0.20 | 1.00 | 1.00 | 0.01 | 0.00 | 1.00 | 1.00 | ok |
| Sc3_rho0 | BART-PP | 0.02 | 0.03 | 0.012 | 0.28 | 0.28 | 0.30 | 0.27 | 0.93 | 0.94 | 0.58 | 0.96 | ok |

## Arm-specific metrics (paper / ours)

| Scenario | Method | Trt bias paper | ours | Trt RMSE paper | ours | Ctrl bias paper | ours | SE | Ctrl RMSE paper | ours |
|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | -0.01 | -0.00 | 0.24 | 0.15 | 0.02 | 0.00 | 0.011 | 0.36 | 0.25 |
| Sc1 | LM-CP | -0.01 | -0.00 | 0.24 | 0.15 | 0.01 | -0.01 | 0.008 | 0.22 | 0.18 |
| Sc1 | LM-PP | -0.01 | -0.00 | 0.24 | 0.15 | 0.00 | 0.00 | 0.011 | 0.36 | 0.24 |
| Sc1 | BART-NP | 0.00 | 0.00 | 0.15 | 0.11 | 0.01 | 0.00 | 0.008 | 0.24 | 0.18 |
| Sc1 | BART-CP | 0.00 | 0.01 | 0.15 | 0.11 | 0.00 | -0.00 | 0.004 | 0.13 | 0.10 |
| Sc1 | BART-PP | 0.00 | 0.01 | 0.15 | 0.11 | 0.00 | -0.00 | 0.008 | 0.23 | 0.17 |
| Sc2 | LM-NP | -0.01 | -0.00 | 0.22 | 0.14 | 0.03 | 0.02 | 0.011 | 0.34 | 0.24 |
| Sc2 | LM-CP | -0.01 | 0.00 | 0.22 | 0.14 | -0.38 | -0.38 | 0.010 | 0.46 | 0.44 |
| Sc2 | LM-PP | -0.01 | 0.00 | 0.22 | 0.14 | -0.03 | -0.01 | 0.011 | 0.36 | 0.25 |
| Sc2 | BART-NP | 0.00 | 0.00 | 0.15 | 0.11 | 0.02 | 0.01 | 0.008 | 0.23 | 0.18 |
| Sc2 | BART-CP | 0.00 | 0.00 | 0.15 | 0.12 | -0.05 | -0.06 | 0.006 | 0.19 | 0.15 |
| Sc2 | BART-PP | 0.00 | 0.00 | 0.15 | 0.12 | 0.00 | -0.00 | 0.007 | 0.23 | 0.16 |
| Sc3_rho0 | LM-NP | 0.00 | 0.02 | 0.25 | 0.15 | 0.03 | 0.01 | 0.011 | 0.38 | 0.26 |
| Sc3_rho0 | LM-CP | 0.00 | 0.02 | 0.25 | 0.15 | -1.01 | -1.02 | 0.008 | 1.02 | 1.04 |
| Sc3_rho0 | LM-PP | 0.00 | 0.02 | 0.25 | 0.15 | -0.06 | -0.03 | 0.011 | 0.37 | 0.24 |
| Sc3_rho0 | BART-NP | 0.02 | 0.02 | 0.20 | 0.13 | 0.03 | 0.02 | 0.010 | 0.33 | 0.22 |
| Sc3_rho0 | BART-CP | 0.02 | 0.02 | 0.20 | 0.13 | -0.96 | -0.96 | 0.006 | 0.96 | 0.97 |
| Sc3_rho0 | BART-PP | 0.02 | 0.02 | 0.20 | 0.13 | 0.00 | -0.01 | 0.009 | 0.30 | 0.20 |

## Largest discrepancies

| Scenario | Method | d bias | d SD | d RMSE | d cov | d pow | d ctrl bias | d ctrl RMSE | d trt RMSE |
|---|---|---|---|---|---|---|---|---|---|
| Sc1 | BART-CP | 0.01 | 0.00 | -0.01 | 0.04 | 0.05 | -0.00 | -0.03 | -0.04 |
| Sc2 | BART-CP | 0.01 | -0.00 | 0.01 | -0.04 | 0.07 | -0.01 | -0.04 | -0.03 |
| Sc3_rho0 | LM-NP | 0.04 | 0.01 | -0.01 | 0.02 | 0.41 | -0.02 | -0.12 | -0.10 |
| Sc1 | LM-CP | 0.03 | 0.00 | -0.02 | 0.03 | 0.37 | -0.02 | -0.04 | -0.09 |
| Sc1 | LM-NP | 0.03 | 0.00 | -0.01 | -0.01 | 0.42 | -0.02 | -0.11 | -0.09 |

## Notes

SE is the Monte Carlo standard error of our bias estimate (R = 500 replicates); the paper's R = 500 values carry an SE of about 0.015 (linear models) to 0.010 (BART). A cell flagged OUT whose discrepancy is within about two combined standard errors is consistent with Monte Carlo variation. All methods within a scenario share the replicate data (common random numbers), so their discrepancies move together and do not count as independent evidence.

## Diagnosis

### Arm-specific RMSE

Our arm-specific RMSE (classical Monte Carlo RMSE of the posterior mean against the truth standardized to the replicate's own RCT profiles) is systematically below the paper's. The table recomputes it under the four candidate definitions: classical or draw-based (the paper's metric-table formula sqrt(B^-1 sum_b (tau_b - tau_true)^2), which adds the posterior variance), against the sample-standardized truth or the population truth E[mu_a(X) | D = 1] (Monte Carlo over 1e5 pool subjects, which adds the covariate-sampling variance of the 300 profiles). Mean absolute deviation from the paper's 36 arm cells:

- classical, sample truth (used above): 0.070
- classical, population truth: 0.034
- draw-based, sample truth: 0.004
- draw-based, population truth: 0.032

| Scenario | Method | Trt paper | cls/smp | cls/pop | drw/smp | drw/pop | Ctrl paper | cls/smp | cls/pop | drw/smp | drw/pop |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | 0.24 | 0.15 | 0.19 | 0.24 | 0.26 | 0.36 | 0.25 | 0.28 | 0.36 | 0.38 |
| Sc1 | LM-CP | 0.24 | 0.15 | 0.19 | 0.24 | 0.26 | 0.22 | 0.18 | 0.16 | 0.22 | 0.21 |
| Sc1 | LM-PP | 0.24 | 0.15 | 0.19 | 0.24 | 0.26 | 0.36 | 0.24 | 0.28 | 0.37 | 0.38 |
| Sc1 | BART-NP | 0.15 | 0.11 | 0.17 | 0.15 | 0.18 | 0.24 | 0.18 | 0.22 | 0.24 | 0.27 |
| Sc1 | BART-CP | 0.15 | 0.11 | 0.17 | 0.15 | 0.18 | 0.13 | 0.10 | 0.16 | 0.13 | 0.17 |
| Sc1 | BART-PP | 0.15 | 0.11 | 0.17 | 0.15 | 0.18 | 0.23 | 0.17 | 0.21 | 0.24 | 0.26 |
| Sc2 | LM-NP | 0.22 | 0.14 | 0.22 | 0.23 | 0.27 | 0.34 | 0.24 | 0.28 | 0.34 | 0.37 |
| Sc2 | LM-CP | 0.22 | 0.14 | 0.22 | 0.22 | 0.27 | 0.46 | 0.44 | 0.46 | 0.46 | 0.48 |
| Sc2 | LM-PP | 0.22 | 0.14 | 0.22 | 0.23 | 0.27 | 0.36 | 0.25 | 0.29 | 0.37 | 0.39 |
| Sc2 | BART-NP | 0.15 | 0.11 | 0.20 | 0.15 | 0.21 | 0.23 | 0.18 | 0.24 | 0.24 | 0.27 |
| Sc2 | BART-CP | 0.15 | 0.12 | 0.20 | 0.15 | 0.20 | 0.19 | 0.15 | 0.22 | 0.20 | 0.24 |
| Sc2 | BART-PP | 0.15 | 0.12 | 0.20 | 0.15 | 0.21 | 0.23 | 0.16 | 0.23 | 0.23 | 0.27 |
| Sc3_rho0 | LM-NP | 0.25 | 0.15 | 0.21 | 0.25 | 0.29 | 0.38 | 0.26 | 0.29 | 0.38 | 0.40 |
| Sc3_rho0 | LM-CP | 0.25 | 0.15 | 0.21 | 0.25 | 0.29 | 1.02 | 1.04 | 1.02 | 1.03 | 1.02 |
| Sc3_rho0 | LM-PP | 0.25 | 0.15 | 0.21 | 0.25 | 0.29 | 0.37 | 0.24 | 0.28 | 0.37 | 0.39 |
| Sc3_rho0 | BART-NP | 0.20 | 0.13 | 0.20 | 0.20 | 0.24 | 0.33 | 0.22 | 0.26 | 0.32 | 0.34 |
| Sc3_rho0 | BART-CP | 0.20 | 0.13 | 0.20 | 0.20 | 0.24 | 0.96 | 0.97 | 0.96 | 0.97 | 0.96 |
| Sc3_rho0 | BART-PP | 0.20 | 0.13 | 0.20 | 0.20 | 0.24 | 0.30 | 0.20 | 0.24 | 0.29 | 0.32 |

### Power

With a true effect of 1 and posterior SDs of 0.14 to 0.37, a 95% credible interval excludes 0 in essentially every replicate, so the rejection rate under the paper's stated definition (Table ate_metrics: 1{0 not in CI95}) and stated effect (delta = 1) is close to 1 for every method. The paper's power column (0.39 to 0.95) cannot arise from that definition at delta = 1. The recomputation below shifts our posterior draws so that the effect is 0.5 (the posterior of a location-shift estimator moves with the effect, so this is a close proxy for a run at delta = 0.5) and applies three criteria. Mean absolute deviation from the paper's 18 power cells:

- 95% CrI excludes 0, effect 1 (paper's definition): 0.259
- 95% CrI excludes 0, draws shifted to effect 0.5: 0.079
- 90% CrI excludes 0, draws shifted to effect 0.5: 0.012
- one-sided 95%, draws shifted to effect 0.5: 0.012

| Scenario | Method | Pow paper | ATE SD ours | 95%, effect 1 | 95%, effect 0.5 | 90%, effect 0.5 | one-sided 95%, effect 0.5 |
|---|---|---|---|---|---|---|---|
| Sc1 | LM-NP | 0.41 | 0.34 | 0.83 | 0.32 | 0.41 | 0.41 |
| Sc1 | LM-CP | 0.62 | 0.24 | 0.99 | 0.53 | 0.65 | 0.65 |
| Sc1 | LM-PP | 0.40 | 0.35 | 0.85 | 0.29 | 0.40 | 0.40 |
| Sc1 | BART-NP | 0.77 | 0.20 | 1.00 | 0.68 | 0.77 | 0.77 |
| Sc1 | BART-CP | 0.95 | 0.14 | 1.00 | 0.96 | 0.98 | 0.98 |
| Sc1 | BART-PP | 0.78 | 0.21 | 1.00 | 0.69 | 0.77 | 0.77 |
| Sc2 | LM-NP | 0.44 | 0.32 | 0.87 | 0.36 | 0.45 | 0.45 |
| Sc2 | LM-CP | 0.91 | 0.29 | 0.99 | 0.84 | 0.90 | 0.90 |
| Sc2 | LM-PP | 0.48 | 0.34 | 0.84 | 0.34 | 0.43 | 0.43 |
| Sc2 | BART-NP | 0.77 | 0.20 | 1.00 | 0.67 | 0.77 | 0.77 |
| Sc2 | BART-CP | 0.93 | 0.18 | 1.00 | 0.87 | 0.92 | 0.92 |
| Sc2 | BART-PP | 0.79 | 0.21 | 1.00 | 0.70 | 0.79 | 0.79 |
| Sc3_rho0 | LM-NP | 0.39 | 0.37 | 0.80 | 0.29 | 0.41 | 0.41 |
| Sc3_rho0 | LM-CP | 1.00 | 0.26 | 1.00 | 1.00 | 1.00 | 1.00 |
| Sc3_rho0 | LM-PP | 0.49 | 0.36 | 0.85 | 0.34 | 0.47 | 0.47 |
| Sc3_rho0 | BART-NP | 0.52 | 0.29 | 0.92 | 0.41 | 0.52 | 0.52 |
| Sc3_rho0 | BART-CP | 1.00 | 0.20 | 1.00 | 1.00 | 1.00 | 1.00 |
| Sc3_rho0 | BART-PP | 0.58 | 0.28 | 0.96 | 0.51 | 0.60 | 0.60 |

### Reading

The ATE bias, posterior SD, RMSE and coverage columns reproduce the paper within Monte Carlo error for all 18 rows once the ATE RMSE is read as the Monte Carlo RMSE of the posterior mean; the data-generating mechanism and the six reference comparators are therefore implemented as the paper describes. Two columns do not reproduce under the paper's stated definitions, and the recomputations above identify what the published numbers are. The arm-specific RMSE matches the draw-based formula of the metric table against the sample-standardized truth (mean absolute deviation 0.004 over 36 cells, versus 0.07 for the classical RMSE used for the ATE), so the paper mixes the two formulas: classical for the ATE, draw-based for the arms. The power column matches a 90% credible interval (equivalently a one-sided 95% test) at an effect of 0.5 (mean absolute deviation 0.012 over 18 cells), and not the stated 95% interval at delta = 1, under which every method has power near 1. We have not tuned anything to close either gap; the harness stores both RMSE variants and the per-replicate draws, so the revision can state its definitions explicitly and recompute power at whatever effect size and level it adopts.
