# ess_local (survival)

# Default (RCT trt n1 = 200, RWD n3 = 300) -- outcome n200, Ntarget 200

## sc1

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 1 --size_option default
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n200_sc1E_alternative_1.RData \
    --outcome n200 --Ntarget 200 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.0001,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n200_nt10_k2_sc1E --mult "1,0.75,0.5"
Rscript ess_plot_simple.R --scenario n200_nt5_k2_sc1E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n200_nt5_k2.png
```

## sc2

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 2 --size_option default
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n200_sc2E_alternative_1.RData \
    --outcome n200 --Ntarget 200 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.0001,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n200_nt5_k2_sc2E --mult "3,2,1,0.75,0.5"
Rscript ess_plot_simple.R --scenario n200_nt50_k2_sc2E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n200_nt50_k2.png
```

## sc3

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 3 --size_option default
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n200_sc3E_cor0.7_alternative_1.RData \
    --outcome n200 --Ntarget 200 --seed 6 --ntree 50 --k 2 \
    --s2-grid "0.0001,0.25" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n200_nt5_k2_sc3E_cor0.7 --mult "0.2,0.1"
Rscript ess_plot_simple.R --scenario n200_nt10_k2_sc3E_cor0.7 --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n200_nt10_k2_sc3.png
```

# Small (RCT trt n1 = 30, RWD n3 = 300) -- outcome n30, Ntarget 30

## sc1

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 1 --size_option small
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n30_sc1E_alternative_1.RData \
    --outcome n30 --Ntarget 30 --seed 6 --ntree 5 --k 2 \
    --s2-grid "0.0001,0.10,0.11,0.15,0.16,0.25,0.26,0.30,0.31,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n30_nt5_k2_sc1E --mult "3,2,1,0.75,0.5"
Rscript ess_plot_simple.R --scenario n30_nt50_k2_sc1E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n30_nt50_k2.png
```

## sc2

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 2 --size_option small
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n30_sc2E_alternative_1.RData \
    --outcome n30 --Ntarget 30 --seed 6 --ntree 5 --k 2 \
    --s2-grid "0.0001,0.15,0.16,0.25,0.26,0.30,0.31,0.38,0.39,0.40,0.42,0.44" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n30_nt5_k2_sc2E --mult "3,2,1,0.75,0.5"
Rscript ess_plot_simple.R --scenario n30_nt50_k2_sc2E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n30_nt50_k2.png
```

## sc3

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/ess_local
Rscript gen_one.R --sc 3 --size_option small
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-survival-single-arm/data_v3/data_p10_n30_sc3E_cor0.7_alternative_1.RData \
    --outcome n30 --Ntarget 30 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.0001,0.15,0.16,0.25,0.26,0.30,0.31,0.38,0.39,0.40,0.42,0.44" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --outcome surv --scenario n30_nt5_k2_sc3E_cor0.7 --mult "3,2,1,0.75,0.5"
Rscript ess_plot_simple.R --scenario n30_nt5_k2_sc3E_cor0.7 --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_surv_real_n30_nt5_k2_sc3.png
```
