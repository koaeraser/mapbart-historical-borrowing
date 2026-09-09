# mapbart-case-study-mm — ESS calibration + analysis pipeline

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm
```

## 1. Clean + merge

```bash
Rscript data_cleaning_ucmm.R     # -> data_cleaned/ucmm_cleaned_n253.RData, _n48.RData
Rscript data_cleaning_elokrd.R   # -> data_cleaned/elokrd_cleaned_n30.RData
Rscript data_merge.R             # -> merged_elokrd_ucmm_n283.RData (primary), _n78.RData (secondary)
```

## 2. Calibrate the prior s² (`ess_local`)

See `ess_local/README.md`.

## 3. Run the analysis pipeline

`MERGED_FILE` may be a bare filename (resolved against `data_cleaned/`) or a full
path. Run `data_merge.R` first so the file exists.

```bash
OUTCOME=PFS MERGED_FILE=merged_elokrd_ucmm_n283.RData MBART_CONFIGS="5:2,10:2,50:2" HIERAFT_CONFIGS="0.05,0.5" Rscript run_all.R
OUTCOME=OS MERGED_FILE=merged_elokrd_ucmm_n283.RData MBART_CONFIGS="5:2,10:2,50:2" HIERAFT_CONFIGS="0.05,0.5" Rscript run_all.R

OUTCOME=PFS MERGED_FILE=merged_elokrd_ucmm_n78.RData MBART_CONFIGS="5:2,10:2,50:2" HIERAFT_CONFIGS="0.05,0.5" Rscript run_all.R
OUTCOME=OS MERGED_FILE=merged_elokrd_ucmm_n78.RData MBART_CONFIGS="5:2,10:2,50:2" HIERAFT_CONFIGS="0.05,0.5" Rscript run_all.R

# Rerun a subset of scripts
SCRIPTS="mBART_analysis.R" Rscript run_all.R

# Rebuild the summary table from existing res/ files (no re-run)
SKIP_RUN=1 Rscript run_all.R
```

Notes:
- `HIERAFT_CONFIGS` is the discrepancy-prior s² sweep (default `0.05,0.5`,
  matching `mapbart-sim-survival-realcov`).
- AFTv2 is run at a sweep of power-prior weights `rwd_w` (external-control
  downweighting). `run_all.R` derives these automatically from the ESS
  calibration: `rwd_w = 1` (full borrowing) plus `N_target / n_UCMM` for each
  MAP-BART target N, so AFTv2 is compared at matched effective control size.
  Override with `AFTV2_RWD_W="1,0.5,..."` if running `AFTv2_analysis.R` directly.
