# Numerical and source audit

## Existing table rows
- tab_app.tex: PASS; every original numeric row retained.
- tab_gauss.tex: PASS; original numerical rows retained. Multirow spans increased for added CAHB rows.
- tab_gauss_full.tex: PASS; every original numeric row retained.
- tab_ladder.tex: PASS; every original numeric row retained.
- tab_map.tex: PASS; original numerical rows retained. Multirow spans increased for added CAHB rows.
- tab_map_surv.tex: PASS; every original numeric row retained.
- tab_single.tex: PASS; every original numeric row retained.
- tab_single_full.tex: PASS; every original numeric row retained.
- tab_surv.tex: PASS; every original numeric row retained.
- tab_surv_full.tex: PASS; every original numeric row retained.

## New prose values
30 context-specific cell/rounding checks passed.
- CAHB CSV: Sc1, CAHB, rmse -> 0.258
- CAHB CSV: Sc1, LRC-BART-100, rmse -> 0.175
- CAHB CSV: Sc1, BART-NP, rmse -> 0.223
- CAHB CSV: Sc4_X5_d2, CAHB, bias -> -0.333
- CAHB CSV: Sc4_X5_d2, CAHB, map_R -> 0.49
- CAHB CSV: Sc4_X5_d2, CAHB, map_Rc -> 0.41
- CAHB CSV: Sc5_d1, CAHB, bias -> -0.437
- CAHB CSV: Sc1, CAHB-lit, sd -> 12.226
- CAHB CSV: Sc1, CAHB-lit, rmse -> 0.950
- CAHB CSV: Sc1, CAHB-p10, sd -> 0.542
- CAHB CSV: Sc1, CAHB-p10, rmse -> 0.262
- CAHB CSV: Sc1, CAHB, sd -> 0.410
- Regional CSV: Sc1, R, ess_real_S -> 32.7
- Regional CSV: Sc1, Rc, ess_real_S -> 31.9
- Regional CSV: Sc4_X5_d1, R, ess_real_S -> 15.1
- Regional CSV: Sc4_X5_d1, Rc, ess_real_S -> 12.7
- Regional CSV: Sc4_X5_d2, R, ess_real_S -> 9.5
- Regional CSV: Sc4_X5_d2, Rc, ess_real_S -> 0.44
- Regional CSV: Sc5_d2, R, ess_real_S -> 0.2
- Regional CSV: Sc5_d2, Rc, ess_real_S -> 0.2
- Paired CSV: Sc4_X5_d2, rmse_R, lrc -> 0.751
- Paired CSV: Sc4_X5_d2, rmse_R, pp -> 0.704
- Paired CSV: Sc4_X5_d2, rmse_R, diff -> 0.047
- Paired CSV: Sc4_X5_d2, rmse_R, se_paired -> 0.005
- Paired CSV: Sc4_X5_d2, rmse_R_far, lrc -> 0.747
- Paired CSV: Sc4_X5_d2, rmse_R_far, pp -> 0.723
- Paired CSV: Sc4_X5_d2, rmse_R_far, diff -> 0.024
- Paired CSV: Sc4_X5_d2, rmse_R_far, se_paired -> 0.005
- Paired CSV: Sc4_X5_d2, rmse_allprof_R, diff -> 0.039
- Paired CSV: Sc4_X5_d2, rmse_allprof_R, se_paired -> 0.006

## Scope
The strict prose checker preserves every numerical token and inline math object against post_technical. Existing table rows were compared with the immutable original. This is not a new independent replication of the original 60-number audit or of the simulations. Source files under 01-code, 02-validation and 04-application were not modified.
The JSM deck was read for scope assessment. Its quoted existing table cells have not changed. Its wording still contains the old one-residual-SD and sample-size-limit statements; no deck rewrite was performed in this manuscript prose pass.
