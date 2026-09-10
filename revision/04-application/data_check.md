# Data check: merged EloKRd + UCMM analysis dataset against Table baseline

Date: 2026-09-09. File: `yunxuan-repo/mapbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n283.RData` (object `merged`, 283 rows, 13 columns: subject_id, source, age, sex, race, ethnicity, high_risk_cyto, transplant_off_protocol, os_months, os_status, pfs_months, pfs_status, trt). No missing values. Times are in months and every script in the repo divides by 12; we do the same. The companion file `merged_elokrd_ucmm_n78.RData` is the secondary cohort of `data_cleaning_ucmm.R` (lines 231 to 240): the same 30 EloKRd patients with UCMM restricted to the KRd and Rd regimens (48 patients, "Cohort 3"). The paper's tables use the n = 283 primary cohort only, and so do we.

## Agreement with Table baseline

| Quantity | Paper | Data | Status |
|---|---|---|---|
| n EloKRd / UCMM | 30 / 253 | 30 / 253 | match |
| Age, median (range), EloKRd | 62.8 (43 to 81) | 62.8 (42.8 to 81.2) | match |
| Age, median (range), UCMM | 61.3 (32 to 97) | 61.3 (31.8 to 97.3) | match |
| Male | 21 (70%) / 130 (51%) | 21 / 130 | match |
| Race White | 23 (77%) / 164 (65%) | 23 / 164 | match |
| Race Black | 4 (13%) / 76 (30%) | 4 / 76 | match |
| Race Other/Unknown | 3 (10%) / 13 (5%) | 3 / 13, but see note 1 | match in count |
| Hispanic or Latino | 1 (3%) / 13 (5%) | 1 / 13 | match |
| Non-Hispanic | 29 (97%) / 240 (95%) | 29 / 240, but see note 2 | match in count |
| High-risk cytogenetics | 15 (50%) / 58 (23%) | 15 / 58 | match |
| ASCT off protocol | 9 (30%) / 203 (80%) | 9 / 203 | match |
| PFS events | 13 (43%) / 159 (63%) | 13 / 159 | match |
| Median PFS, years | not reached / 5.08 | not reached / 5.087 | 5.087 rounds to 5.09, not 5.08 |
| PFS RMST(5 yr), years | 3.71 / 3.65 | 3.705 / 3.649 | match |
| OS events | 10 (33%) / 60 (24%) | 10 / 60 | match |
| Median OS | not reached / not reached | not reached / not reached | match |
| OS RMST(5 yr), years | 4.12 / 4.57 | 4.122 / 4.572 | match |
| Log-rank p, PFS / OS | 0.506 / 0.048 | 0.506 / 0.0476 | match |
| KM 5-year RMST ratio, PFS | 1.015 | 3.705 / 3.649 = 1.015 | match |

Maximum follow-up is 6.75 years in EloKRd and 24.1 years in UCMM, so the 5-year horizon lies inside both, as the paper states.

## Mismatches and coding issues

1. Race has four levels in the UCMM rows (Asian 4, Black 76, Other/Unknown 9, White 164) and three in EloKRd. Table baseline reports 13 Other/Unknown for UCMM, so it pools Asian with Other/Unknown, but the analysis scripts (`prep_covariates()` in every `*_analysis.R`) build dummies from the raw four-level factor, giving race_Black, race_Other_Unknown and race_White with Asian as the reference level. The count in the table is right; the covariate in the models is not the covariate in the table.

2. Ethnicity is not harmonized. The non-Hispanic label is "Non-Hispanic" in the 29 EloKRd rows and "Not Hispanic or Latino" in the 240 UCMM rows. `prep_covariates()` therefore makes two dummies: `ethnicity_Non-Hispanic` (1 for 29 EloKRd patients, 0 for all UCMM) and `ethnicity_Not_Hispanic_or_Latino` (1 for 240 UCMM patients, 0 for all EloKRd). In every model-based row of Table rmst-results the control model is trained on UCMM alone and evaluated at the 30 EloKRd profiles, so (a) the first dummy is constant zero in training (dropped by `bartModelMatrix(rm.const = TRUE)` for the BART rows; retained with a prior-only coefficient in the Stan AFT rows, where it adds a N(0, lambda_beta^2) draw to the predicted log time of the 29 non-Hispanic EloKRd patients), and (b) the 29 non-Hispanic EloKRd patients carry `ethnicity_Not_Hispanic_or_Latino = 0`, the same code as the 13 Hispanic UCMM patients, so the counterfactual control prediction for them is the UCMM Hispanic-stratum prediction. Table baseline hides this because it reports the counts, which agree. The paper's reference rows and its MAP-AFT-BART rows were all run under this coding.

3. Minor: the UCMM median PFS is 5.087 years, printed as 5.08 in the table (5.09 at two decimals).

## What we use downstream

We run the reference models under both codings: the original `prep_covariates()` coding, to reproduce the published rows, and a harmonized seven-column coding (age, male, race_Black, race_Other with Asian pooled into Other/Unknown as in Table baseline, hispanic, high_risk_cyto, asct), which is what the six covariates of the paper describe. LRC-BART is fitted under the harmonized coding, with the original coding as a sensitivity so that its effect can be separated from the change of method. Code: `revision/04-application/code/prep_data.R`.
