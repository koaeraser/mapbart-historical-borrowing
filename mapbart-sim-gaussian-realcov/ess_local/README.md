# ess_local (Gaussian)

`ess_plot.R` and `ess_plot_simple.R` append the scenario tag
(`sc1E` / `sc2E`, parsed from the calibration's source RWD file) to the
output filename, so figures from different scenarios never overwrite each
other. Run each scenario one by one: calibration, then the two figure steps,
before moving on to the next.

Note: `--ntree`/`--k` set the control-arm BART config; the calibration `.RData`,
the plot `--scenario`, and the `--out` PNG are all tagged `_nt<ntree>_k<k>` so
configs never clobber each other. Defaults are `--ntree 50 --k 2`. Calibrate
each config the MAP-BART sweep will use.

# Default RWD (n3 = 300) -- Ntarget 200

## sc1

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian-realcov/ess_local
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian-realcov/data_v3/data_p10_sc1E_alternative_1.RData \
    --Ntarget 200 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --scenario nt10_k2_sc1E --mult "1,0.75,0.5"
Rscript ess_plot_simple.R --scenario nt5_k2_sc1E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_gaussian_real_nt5_k2.png
```

## sc2

```bash
cd /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian-realcov/ess_local
CORES=1 Rscript ess_calibration.R \
    --data /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian-realcov/data_v3/data_p10_sc2E_alternative_1.RData \
    --Ntarget 200 --seed 6 --ntree 10 --k 2 \
    --s2-grid "0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11" \
    --map-model synthetic --beta-prior auto --sigma-ref auto --B 100 --burn 1000 --q 0.95
Rscript ess_plot.R --scenario nt10_k2_sc2E --mult "1,0.75,0.5"
Rscript ess_plot_simple.R --scenario nt50_k2_sc2E --out /Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-sim-gaussian/manuscript/inserts/ess_gaussian_real_nt50_k2.png
```

To calibrate a different control config, change `--ntree`/`--k` (e.g. `--ntree 10
--k 2`); every output is tagged accordingly (`..._nt10_k2_sc<...>`). Run one config
per invocation, and calibrate every config the MAP-BART sweep will use.
