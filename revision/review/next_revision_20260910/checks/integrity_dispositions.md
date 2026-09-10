# Integrity audit dispositions

This run includes verified scientific changes, so the original-to-final humanize-prose checker is not expected to be a prose-only PASS. Its hard failures were inspected rather than suppressed. The separate post-technical-to-shortened-main checks in round_1/prose_s12.log and prose_s35.log passed hard preservation checks; their inline-math removal audit concerns repeated symbols, not changed definitions.

| Original-to-final hard delta | Disposition |
|---|---|
| sections_1_2 equations | Equation (6) terminal comma changed to a period; its numerical calibration formula is unchanged. Unnumbered variance formula and inline theorem conditions are documented scientific corrections. |
| sections_1_2 references | Added precise calibration/Pr-rule cross-references and removed repeated model reference during shortening. Targets remain defined. |
| sections_3_5 references | Removed a repeated single-arm cross-reference from the replaced false unchanged-survival-mean argument. |
| web_appendix_G sections | G.6 title now states conditional prior-cost scope. |
| web_appendix_G references | Redundant forward/backward pointers removed with replaced model/coding/tail explanations. The referenced sections/tables remain present. |
| appendix_new labels | Added eqA:leafvariance, defining the exact leaf-weighted variance. No existing label removed. |
| appendix_new sections | D.1 added; D renamed to defined variance-based calibration; C and F headings replaced development history with scientific scope. |
| appendix_new equations | Added leaf variance; ESS/Pr-rule punctuation; positivity condition expressed as a complete equivalence. Appendix D's unnumbered theorem/proof equations replaced as recorded in both review ledgers. |
| appendix_new references | Removed obsolete proof/report/old-theorem accounting; tied current calibration, leaf variance and simulations to their correct labels. |
| slides citations | Added the verified Jin2023CAHB citation already available in the bibliography. |

Number/number-word, inline-math and structure audits were inspected in the associated source changes. Material changes are the corrected leaf-state sum, working-ESS scope, explicit cap/fallback, rescaled block conditions, finite-prior anchor, conditional G.6 optimizer, nonlinear survival examples, correct application cohort/ASCT and sigma assumptions, and removal of unsupported 17–86 threshold/obsolete 0.7 illustration. Existing numerical result corrections are traced to 19 producing CSV rows in quoted_corrections.json. Rhetorical repetitions were shortened without dropping scenarios or results. Main theorems now require interior mixing weight. G.6's false optimizer was replaced, not preserved as an invariant.

The final side-review integration replaces empirical “unbiased” wording and clarifies counterfactual variance equality and finite-grid selection. The user's indicator question prompted an explicit all-spike label and z-to-variance explanation in slide 6. Main wrapper change only hides hyperlink borders. Slide image sizing/reference line spacing changes prevent overflow while retaining the 25-slide structure. No numerical table or figure files changed. No undefined references, citations, duplicate-label warnings or overfull boxes remain in the manuscript builds.
