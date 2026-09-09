# ess_local (mapbart-case-study-mm)

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm/ess_local
```

Note: `--ntree`/`--k` set the AFT-BART control config; the calibration `.RData`
(and its checkpoint) are tagged `_nt<ntree>_k<k>` so configs never clobber each
other, and `mBART_analysis.R` matches the calibration by this tag. Defaults are
`--ntree 50 --k 2`. Run each cohort/outcome/config one by one: calibration, then
the figure step.

## 1. Calibration

```bash
# merged_elokrd_ucmm_n283.RData (primary: all regimens, E-Rd excluded)
# PFS
CORES=1 Rscript ess_calibration.R --data merged_elokrd_ucmm_n283.RData \
    --outcome PFS --Ntarget 30 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.08,0.09,0.10,0.11,0.12,0.15,0.16,0.28,0.29" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95

# OS
CORES=1 Rscript ess_calibration.R --data merged_elokrd_ucmm_n283.RData \
    --outcome OS  --Ntarget 30 --seed 6 --ntree 50 --k 2 \
    --s2-grid "0.28,0.29,0.3,0.33,0.34,0.35,0.42,0.43,0.44" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95

# merged_elokrd_ucmm_n78.RData (secondary: KRd/Rd only)
# PFS
CORES=1 Rscript ess_calibration.R --data merged_elokrd_ucmm_n78.RData \
    --outcome PFS --Ntarget 30 --seed 6 --ntree 50 --k 2 \
    --s2-grid "0.25,0.26,0.27,0.29,0.30,0.31,0.38,0.39,0.40" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95

# OS
CORES=1 Rscript ess_calibration.R --data merged_elokrd_ucmm_n78.RData \
    --outcome OS  --Ntarget 30 --seed 6 --ntree 50 --k 2 \
    --s2-grid "0.28,0.29,0.3,0.33,0.34,0.35,0.42,0.43,0.44" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
```

## 2. Figures

The `--res` path must include the config tag written by the calibration.

```bash
Rscript ess_plot.R --outcome surv --mult "3,2,1,0.75,0.5" --res res/ess_calibration_ucmm_pfs_synthetic_n283_nt10_k2.RData
Rscript ess_plot.R --outcome surv --mult "3,2,1,0.75,0.5" --res res/ess_calibration_ucmm_os_synthetic_n283_nt50_k2.RData
Rscript ess_plot.R --outcome surv --mult "3,2,1,0.75,0.5" --res res/ess_calibration_ucmm_pfs_synthetic_n78_nt50_k2.RData
Rscript ess_plot.R --outcome surv --mult "3,2,1,0.75,0.5" --res res/ess_calibration_ucmm_os_synthetic_n78_nt50_k2.RData
```
