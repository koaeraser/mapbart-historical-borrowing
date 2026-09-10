| # | scenario | criterion | value | result |
|---|---|---|---|---|
| 1 | Sc3_rho-0.5 | |bias| <= 0.05 | 0.064 | FAIL |
| 2 | Sc3_rho-0.5 | coverage >= 0.92 | 0.930 | PASS |
| 3 | Sc3_rho-0.5 | RMSE <= BART-PP + 0.02 | 0.270 vs 0.261 + 0.02 | PASS |
| 4 | Sc3_rho0 | |bias| <= 0.05 | 0.062 | FAIL |
| 5 | Sc3_rho0 | coverage >= 0.92 | 0.940 | PASS |
| 6 | Sc3_rho0 | RMSE <= BART-PP + 0.02 | 0.269 vs 0.262 + 0.02 | PASS |
| 7 | Sc3_rho0.5 | |bias| <= 0.05 | 0.063 | FAIL |
| 8 | Sc3_rho0.5 | coverage >= 0.92 | 0.955 | PASS |
| 9 | Sc3_rho0.5 | RMSE <= BART-PP + 0.02 | 0.266 vs 0.254 + 0.02 | PASS |
| 10 | Sc1 | RMSE <= 0.90 x BART-PP | 0.170 vs 0.90 x 0.206 = 0.185 | PASS |
| 11 | Sc1 | power >= BART-PP + 0.05 | 1.000 vs 1.000 + 0.05 | FAIL |
| 12 | Sc4_X5_d2 | |bias| <= 0.10 | -0.026 | PASS |
| 13 | Sc4_X5_d2 | control RMSE in R <= 0.90 x BART-PP's | 0.815 vs 0.90 x 0.780 = 0.702 | FAIL |
| 14 | Sc4_X5_d2 | all-spike map mean in R >= 0.5 | 0.296 | FAIL |
| 15 | Sc4_X5_d2 | all-spike map mean in R^c <= 0.15 | 0.023 | PASS |
| 16 | Sc4_X5_d2 | P(|g| < 0.5) mean in R >= 0.7 | 0.654 | FAIL |
| 17 | Sc4_X5_d2 | P(|g| < 0.5) mean in R^c <= 0.15 | 0.057 | PASS |
| 18 | Sc4_X5_d2 | mean g within 0.25 of 0 in R | -0.262 | FAIL |
| 19 | Sc4_X5_d2 | mean g within 0.25 of -2 in R^c | -1.682 | FAIL |
| 20 | Sc4_X5_d1 | |bias| <= 0.10 | -0.071 | PASS |
| 21 | Sc4_X5_d1 | overall RMSE no worse than BART-PP | 0.230 vs 0.203 | FAIL |
| 22 | Sc5_d2 | |bias| <= 0.10 | -0.027 | PASS |
| 23 | Sc5_d2 | RMSE <= 1.10 x BART-NP | 0.203 vs 1.10 x 0.223 = 0.246 | PASS |
| 24 | Sc5_d2 | all-spike map mean <= 0.15 (R empty: R^c mean) | 0.004 | PASS |
| 25 | Sc5_d2 | P(|g| < 0.5) mean <= 0.15 | 0.009 | PASS |

PASS 15 of 25
