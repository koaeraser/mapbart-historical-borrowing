# Simulation infrastructure for the LRC-BART revision (Phase 1a)

This directory holds the data-generating mechanisms, the paper's comparator methods, the metrics, and a parallel resumable harness for the Biometrics revision of the MAP-BART paper. The proposed method (LRC-BART) is developed separately and plugs in through the method registry (see "Adding a method").

## Layout

- `R/dgp.R`: `gen_data()` for Sc1 to Sc5 (Gaussian and survival), `parse_scenario()` and `gen_data_id()` for the string scenario ids, the survival helpers `lnorm_rmst()` and `lnorm_pop_median()`.
- `R/estimands.R`: standardization of posterior draws of the arm surfaces at the RCT profiles into the paper's estimands (ATE; population median ratio and RMST ratio at tau = 3 for survival), and `make_fit()`, the common fit object.
- `R/comparators.R`: `fit_<name>(dat, outcome, ...)` for LM-NP/CP/PP, BART-NP/CP/PP, PSCL, MAP (N_target 100/75/50), HierLM, and the survival analogues AFT-NP/CP/PP, HierAFT; the registry `METHODS`, `register_method()`, `default_methods()`. The header comment lists every deviation from the paper.
- `R/metrics.R`: per-replicate metrics (`compute_metrics()`), aggregation (`aggregate_metrics()`), and the Table 1 view.
- `R/harness.R`: `run_sim()`, `summarize_sim()`, `time_methods()`, `source_all()`.
- `stan/hier_lm.stan`: HMC version of HierLM used only to cross-check the Gibbs sampler.
- `tests/test_dgp.R`: DGP checks against the paper's summary facts. `tests/check_hier_stan.R`: Gibbs versus Stan cross-check for HierLM.
- `run_reference.R`: reproduces the paper's Table 1 rows for the six reference comparators and writes `reports/reference_check.md`.
- `results/`: simulation output (one RDS per scenario, method, replicate; summaries as CSV). `results/timings.csv` holds the measured per-fit timings.

## Requirements

R 4.4 with packages BART (2.9.10 used), RBesT (1.8.2), survival, data.table, parallel (base), here, testthat, and rstan (only for `tests/check_hier_stan.R`). Parallelism uses forked workers through `parallel::mclapply`, so the harness runs on macOS and Linux; on Windows set `workers = 1`.

## Running

All scripts locate the code root through `here::here()` (a `.here` marker sits in this directory) or a `CODE_ROOT` variable set before sourcing, so run them from anywhere with the working directory anywhere; the examples assume the working directory is this directory.

```
Rscript tests/test_dgp.R            # DGP checks, about one minute
Rscript tests/check_hier_stan.R     # Gibbs versus Stan, a few minutes (compiles Stan)
Rscript run_reference.R 100 18      # reference check, R = 100 replicates on 18 workers
Rscript run_reference.R 500 18      # resumes and extends the same run to R = 500
```

A study is run from R as

```r
CODE_ROOT <- here::here(); source(file.path(CODE_ROOT, "R", "harness.R")); source_all()
run_sim(scenarios = c("Sc1", "Sc2", "Sc3_rho-0.5", "Sc3_rho0", "Sc3_rho0.5",
                      "Sc4_X5_d0.5", "Sc4_X5_d1", "Sc4_X5_d2", "Sc4_X7_d1", "Sc5_d1", "Sc5_d2"),
        outcome = "gaussian", methods = default_methods("gaussian"),
        R = 200, seed = 2026, workers = 18, out_dir = "results/dev_gaussian")
res <- summarize_sim("results/dev_gaussian")   # writes summary_full.csv, summary_table1.csv, per_replicate.csv
res$table1
```

Scenario ids: `Sc1`, `Sc2`, `Sc3_rho<r>` with r in -0.5, 0, 0.5; `Sc4_<region>_d<delta>` with region `X5` (primary, R = {X5 > 2}) or `X7` (secondary, R = {X7 <= 2}) and delta the RWD shift outside R; `Sc5_d<delta>` (shift everywhere); append `_slope` for the optional slope variant. The replicate data for replicate r have seed `seed + r` in every scenario, so the same replicate is paired across methods and scenarios (Sc4 with delta 0 reproduces Sc1 replicate by replicate). Each fit sets its own seed from (seed, r, method name). Existing RDS files are skipped, so an interrupted run resumes by repeating the call; fit errors are written as `.err` files next to the missing RDS and logged in `run.log`.

Per-method timing: `time_methods("Sc2", "gaussian", default_methods("gaussian"), n_rep = 3)`.

## Adding a method

Write `fit_lrcbart <- function(dat, outcome, ...)` that returns `make_fit(mu1, mu0, outcome, sig1, sig0, tau, map = , ess = )`, where `mu1` and `mu0` are B x n_rct matrices of posterior draws of the arm surfaces at `dat$x_rct` (treated profiles first, then control; `dat$z_rct` gives the arm), `sig1`, `sig0` are the residual-scale draws (survival only), `map` is P(borrow) at the RCT control profiles (length `sum(dat$z_rct == 0)`), and `ess` is `list(prior = , realized = )`. Then `register_method("LRC-BART-100", fit_lrcbart)` and pass the name to `run_sim()`. The metrics module picks up `map` and `ess` and produces the Sc4 region-restricted surface RMSE and the borrowing-map summaries automatically.

## Metric definitions

Per replicate and estimand: posterior mean, posterior SD, the draw-based RMSE of the paper's metric table, the central 95% credible interval, coverage, and rejection of the null (0 for the Gaussian ATE, 1 for the survival ratios). Aggregated: bias is the mean of (posterior mean minus truth); `rmse` is the Monte Carlo RMSE of the posterior mean, which is what the paper's tables contain (the draw-based version is also reported as `rmse_draws`; see `R/metrics.R` and `reports/reference_check.md`); coverage and power are proportions. Truths are standardized to the replicate's RCT profiles (`dat$truth`); the population truths are in `dat$truth_pop`. Survival estimands are reported twice, with suffix `_med` (median ratio, the paper's primary) and `_rmst` (RMST ratio at tau = 3).

## Measured timings

Seconds per fit (mean of 3 replicates on Sc2, one core, from `results/timings.csv`; M-series Mac, R 4.4.1). A fit includes the treatment-arm model and the control model, 1000 burn-in and 1000 retained draws.

| Method | Gaussian | Survival |
|---|---|---|
| LM-NP / AFT-NP | 0.07 | 0.25 |
| LM-CP / AFT-CP | 0.07 | 0.25 |
| LM-PP / AFT-PP | 0.07 | 0.26 |
| BART-NP | 0.65 | 0.93 |
| BART-CP | 0.94 | 1.30 |
| BART-PP | 0.96 | 1.28 |
| PSCL | 0.03 | (Gaussian only) |
| MAP-100 / 75 / 50 | 0.55 / 0.61 / 0.63 | (Gaussian only) |
| HierLM / HierAFT | 0.10 | 0.32 |

The full Gaussian comparator set is about 5.5 s per replicate, so 11 scenarios at R = 500 take about 30 core-hours over all comparators, under two hours on 18 workers, before adding the proposed method.

## Reference check

`reports/reference_check.md` compares the six reference comparators with the paper's Table 1 under Sc1, Sc2 and Sc3 (rho = 0) at R = 500 (5 minutes on 18 workers). The ATE bias, SD, RMSE and coverage reproduce within Monte Carlo error (3 of 54 cells outside the 0.03 tolerance, all within two SEs). Two published columns follow other definitions than the paper states: the arm-specific RMSE is the draw-based formula (matched to 0.004), and the power column corresponds to a 90% interval at an effect of 0.5 (matched to 0.012), not the 95% interval at delta = 1, under which power is near 1 for every method. The harness keeps both RMSE variants (`rmse`, `rmse_draws`) and the draws, so the revision can choose and state its definitions.
