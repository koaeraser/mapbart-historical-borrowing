# ess_local (survival)

`ess_plot.R` and `ess_plot_simple.R` append the scenario tag
(`sc1E` / `sc2E` / `sc3E_cor<rho>`, parsed from the calibration's source RWD
file) to the output filename, so figures from different scenarios never
overwrite each other. The calibration `res/` file name is shared across
scenarios, so run each scenario one by one: calibration, then the figure
steps, before moving on to the next.

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival/ess_local
DATA=/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival/data_v2
INSERTS=/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts
```

## sc1 (random RWD selection)

```bash
CORES=8 Rscript ess_calibration.R \
    --data $DATA/data_p10_sc1E_alternative_1.RData \
    --Ntarget 100 --seed 6 \
    --s2-grid "0.18,0.19,0.20,0.21,0.22,0.23,0.28,0.29,0.30" \
    --map-model hybrid --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario sc1E --mult "3,2,1,0.75,0.5"   # -> ..._sc1E.png
Rscript ess_plot_simple.R --scenario sc1E --out $INSERTS/ess_surv.png       # -> ..._sc1E.png
```

## sc2 (measured-confounder RWD selection)

```bash
CORES=8 Rscript ess_calibration.R \
    --data $DATA/data_p10_sc2E_alternative_1.RData \
    --Ntarget 100 --seed 6 \
    --s2-grid "0.18,0.19,0.20,0.21,0.22,0.23,0.28,0.29,0.30" \
    --map-model hybrid --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario sc2E --mult "3,2,1,0.75,0.5"   # -> ..._sc2E.png
Rscript ess_plot_simple.R --scenario sc2E --out $INSERTS/ess_surv.png       # -> ..._sc2E.png
```

## sc3 (unmeasured confounder U) -- one run per correlation rho

### rho = -0.5

```bash
CORES=8 Rscript ess_calibration.R \
    --data $DATA/data_p10_sc3E_cor-0.5_alternative_1.RData \
    --Ntarget 100 --seed 6 \
    --s2-grid "0.18,0.19,0.20,0.21,0.22,0.23,0.28,0.29,0.30" \
    --map-model hybrid --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario sc3E_cor-0.5 --mult "3,2,1,0.75,0.5"   # -> ..._sc3E_cor-0.5.png
Rscript ess_plot_simple.R --scenario sc3E_cor-0.5 --out $INSERTS/ess_surv.png       # -> ..._sc3E_cor-0.5.png
```

### rho = 0

```bash
CORES=8 Rscript ess_calibration.R \
    --data $DATA/data_p10_sc3E_cor0_alternative_1.RData \
    --Ntarget 100 --seed 6 \
    --s2-grid "0.18,0.19,0.20,0.21,0.22,0.23,0.28,0.29,0.30" \
    --map-model hybrid --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario sc3E_cor0 --mult "3,2,1,0.75,0.5"   # -> ..._sc3E_cor0.png
Rscript ess_plot_simple.R --scenario sc3E_cor0 --out $INSERTS/ess_surv.png       # -> ..._sc3E_cor0.png
```

### rho = 0.5

```bash
CORES=8 Rscript ess_calibration.R \
    --data $DATA/data_p10_sc3E_cor0.5_alternative_1.RData \
    --Ntarget 100 --seed 6 \
    --s2-grid "0.18,0.19,0.20,0.21,0.22,0.23,0.28,0.29,0.30" \
    --map-model hybrid --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario sc3E_cor0.5 --mult "3,2,1,0.75,0.5"   # -> ..._sc3E_cor0.5.png
Rscript ess_plot_simple.R --scenario sc3E_cor0.5 --out $INSERTS/ess_surv.png       # -> ..._sc3E_cor0.5.png
```
