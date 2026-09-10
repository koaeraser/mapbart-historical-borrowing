PASS-WITH-FIXES

Scope: sections_3_5.tex, web_appendix_G.tex, tables/*.tex checked against 02-validation/full/*/summary_table.csv, RESULTS_SUMMARY.md, application_report.md, results_table.csv, fig_sc4_map_data.csv and the reconciliation report; cross-references into sections_1_2.tex and 03-theory/appendix_new.tex. Sixty-odd quoted numbers were spot-checked (all of the Sc1-Sc5 Gaussian paragraphs, the survival paragraph, the single-arm paragraph, the application results, G.2-G.5). All match to rounding except the timing claim (Challenge 1); the remaining discrepancies are rounding or wording and are listed below. Every table and figure cited exists and is numbered as cited (Tables 1-5, Figures 1-2, Web Tables S1-S6, Web Figure S1); every comparator named in the text appears in a table; power, RMSE and point-estimate definitions in the text agree with the captions.

## Challenge 1: Timing numbers have no source and contradict the result files
Severity: major
Location: sections_3_5.tex:35
Issue: "A full LRC-BART fit, calibration included, takes 1.2 seconds for Gaussian and 2.1 seconds for survival outcomes on one core of a laptop, against 1.0 and 1.3 seconds for BART-PP." No such benchmark exists in 02-validation or quoted_numbers.txt.
Evidence: the only timing data are the `seconds` columns of summary_table.csv (16-worker run): LRC-BART-100 averages 2.27 s (Gaussian) and 2.54 s (survival); BART-PP 1.62 s and 2.14 s; BART-NP 1.11 s and 1.61 s. The quoted BART-PP times look like BART-NP's, and the LRC-BART times are about half of what was recorded.
Suggested fix: quote the CSV means with their provenance ("per fit, averaged over scenarios, 16 workers"), or run and file a one-core benchmark and cite it.

## Challenge 2: Printed prior ESS exceeds the ceiling, contradicting Theorem 3(a)
Severity: major
Location: tables/tab_single.tex (Table 5), tables/tab_ladder.tex (Web Table S3, Sc2), tables/tab_app.tex (original-coding row), sections_3_5.tex:70
Issue: Section 2.4 and Theorem 3(a) state ESS_0 -> sigma_1^2/V_mu^f as s_0^2 -> 0, so ESS_0 can never exceed the ceiling. The tables print it above the ceiling in every capped row, and a referee will ask why.
Evidence: single_arm/summary_table.csv: LRC-BART (pooled) prior ESS 172.3 with ceiling 163.7 (Gaussian Sc1), 80.2 vs 75.1 (survival n_T=30 Sc1), 34.8 vs 32.8 (Gaussian Sc2); Web Table S3 Sc2 prior ESS 24.5 against ceiling 23; Table 4 original coding ESS_0 12.9 / 15.0 against ceilings 12.6 / 14.6 quoted in the text (line 84). Cause: lrcbart.R:413 sets the ceiling to sigma_1^2/V_whole while line 458 sets ess0_at_s0 to the block mean of sigma_1^2/(V_b + ...), and E[1/V] > 1/E[V].
Suggested fix: report the ceiling as the block mean of sigma_1^2/V^(b) (the quantity the Pr-rule compares against), or cap the printed ESS_0 at the ceiling and add one sentence in Section 2.4 saying the block average is Jensen-biased upward by about 5%.

## Challenge 3: Section 2 promises stronger than Section 3 delivers
Severity: major
Location: sections_1_2.tex:172, 189, 140; sections_3_5.tex:44, 68
Issue: three statements in Sections 2.4-2.6 are not matched by the study. (i) Line 172: "abstention under a clinically sized shift is near certain." Under Sc4 at delta_2 = 1 (0.67 residual SD, a clinically sized shift) the map outside R is 0.52, the ATE bias is -0.041 and coverage 0.89 (Table 2, tab_map.tex); abstention is partial. (ii) Line 189: in single-arm mode "the target is stated as the trial size or a fraction of it"; Section 3.4 uses the absolute target 100 (not n_T = 200 or 30) and fractions of the ceiling, not of the trial size. (iii) Line 140: "In the simulation study c is one third of the RCT residual standard deviation"; that holds for Gaussian outcomes (0.5/1.5) but the survival study uses c = 0.5 with sigma_1 = 1.4 (Web Table S4 caption), 0.36.
Suggested fix: (i) "abstention under a shift of several spike scales is near certain; Section 3 shows what happens at a shift near d_1/2"; (ii) "the target is stated as an absolute number or a fraction of the ceiling"; (iii) "one third of sigma_1 for Gaussian outcomes and 0.5 on the log-time scale for survival outcomes".

## Challenge 4: Sign of the Sc4/Sc5 shift is stated inconsistently
Severity: minor
Location: sections_3_5.tex:29, 51 (Figure 1 caption) versus web_appendix_G.tex:14 and tables/tab_gauss.tex caption
Issue: Section 3.1 and the Figure 1 caption say the RWD surface is shifted "by -delta_2"; G.1 says the RWD mean is f_0 + delta_2 outside R so that g = -delta_2; Table 1's caption says "by delta_2".
Evidence: the negative BART-CP bias under Sc4 (-0.338 at delta_2 = 1, -0.685 at 2; Table 1) requires the RWD control mean to be shifted upward, so G.1 is right and "-delta_2" in the text describes g, not the RWD shift.
Suggested fix: write "the RWD surface is raised by delta_2 outside R, so the true discrepancy g is -delta_2 there" in Section 3.1 and the caption.

## Challenge 5: Discussion limitations omit one of the three stated theory gaps
Severity: minor
Location: sections_3_5.tex:120; appendix_new.tex C.3 "Gap" paragraph
Issue: the appendix lists three unproved items: representability of the discrepancy on the additive span of the g partitions (interactions across two trees are not representable, and the cell means converge to a projection), the posterior over partitions, and f fixed. The Discussion states the second and third and not the first. Web Appendix G.3 also says "five Sc4 delta_2 = 2 criteria" but enumerates four (surface RMSE in R, all-spike map, g in each region); the fifth is P(|g| < 0.5) in R = 0.643 (RESULTS_SUMMARY.md).
Suggested fix: add one clause on representability to the Discussion; list the fifth criterion in G.3.

## Challenge 6: "LRC-BART (100) coincides with the pooled model in every design" overstates
Severity: minor
Location: sections_3_5.tex:70
Issue: the point estimates coincide, the posteriors do not where the target is below the ceiling.
Evidence: Web Table S6, Gaussian Sc1 n_T = 200: LRC-BART (100) SD 0.252 and coverage 1.00 against 0.165 and 0.94 for the pooled model (prior ESS 81.9 vs 172.3).
Suggested fix: "gives the same point estimate as the pooled model in every design (within 0.005) and the same posterior wherever the target is capped".

## Minor observations
- Line 82: unadjusted PFS RMST ratio 1.02 in the text, 1.01 in Table 4 (CSV 1.015); pick one.
- Line 44: Sc5 delta_2 = 2 RMSE "equal to BART-PP's": 0.203 vs 0.205.
- Lines 46 and 118: "abstains at a bias of at most 0.04": Sc3 biases are 0.041-0.042; say "about 0.04".
- Line 42: "one and a half Monte Carlo standard errors": 0.012/0.009 = 1.3.
- Line 44: at delta_2 = 0.5 "g is -0.18 everywhere": Table 2 gives -0.15 in R and -0.18 outside.
- Line 44: "from X_5 = 1.5, where g averages -1.66, to X_5 = 2.5, where it averages -0.44" quotes bin means of (1.5, 2] and (2, 2.5]; say "in the half-unit bins on either side of the boundary".
- Line 98: "H_f = 50 changes the means by at most 0.01": application_report.md says 0.014.
- Line 61: "about 28% censoring": the DGP gives 25% for RCT controls, 29% treated, 27% RWD; "about a quarter" or state which arm.
- web_appendix_G.tex:3 comment says tables S1-S5; there are six (S6 is cited in Table 5's caption).
- figures/fig_rmst_posterior.pdf is built but never cited; either cite it in G.5 or drop it.
- Length: Sections 1-5 run about 10,500 words against a Biometrics norm near 6,500. Candidates to move to Web Appendix G with a one-sentence pointer: the comparator paragraph (line 33, about 330 words, largely duplicated in G.2); the earlier-version reconciliation sentences at lines 61, 70 and 98 (G.2 and G.5 already carry them); the second half of line 84 (coding mechanics); the development-run sentences in line 44 (already in G.3); the prior-cost calculation of Section 2.2 (sections_1_2.tex:82). The first paragraph of the Discussion repeats the Introduction and can go.
- Line 44 is a 600-word paragraph; split at "The control surface is where it costs" and at "At delta_2 = 1".
