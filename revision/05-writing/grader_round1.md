Weighted score: 34/50, grade B

# Grader report, round 1: LRC-BART manuscript (main_draft.tex, 2026-09-09)

We graded the assembled manuscript (sections_1_2.tex, sections_3_5.tex, web_appendix_G.tex, ../03-theory/appendix_new.tex, tables/) against the code in 01-code, the R = 500 results in 02-validation/full, the application in 04-application, and the three earlier reports (verification_phase1.md, go_nogo_report_v2.md, critic_report.md). A Code Auditor subagent checked text against code before we scored Correctness; its table is below.

## Dimension scores

**Correctness: 4/5 (weight 1).** The sampler was verified independently against closed forms, brute-force integration and a two-leaf exact posterior, and the theory critic found zero major and three minor items, all of which the appendix now words correctly (truncated prior stated, A4(a) limit under truncation, aperiodicity). The audit found every hyperparameter the harness ran with matches Section 2 and 3.1. Three text-versus-code wording mismatches remain (Sc4 shift sign in Section 3.1, the surface-RMSE profile set, the "one third of sigma_1" description of c = 0.5 for survival), and the build has seven multiply-defined equation labels because appendix_new.tex reuses eq:model, eq:Mg and others, so cross-references in the compiled PDF may point to the wrong equation.

**Completeness: 3/5 (weight 1).** Section 1.1 motivates Sc4 and Sc5 as "the scenarios a covariate-adaptive method exists for", yet no covariate-adaptive comparator (Jin et al. 2023 CAHB, LEAP, SAM-HC) is run; PSCL is stratified on a propensity score and is not local on the outcome. The power column is uninformative at delta = 1 (1.00 for every tree model) and the design was not changed to a delta at which it discriminates. Sc4 reports no pointwise ESS in R, which is the quantity that would show borrowing where the sources agree, and the realized global ESS of 0.8 at delta_2 = 2 reads as "no borrowing" without it. Paired Monte Carlo standard errors for the LRC-BART versus BART-PP contrasts are not given.

**Rigor: 4/5 (weight 1).** Replicates are paired across methods and scenarios, 25 criteria were pre-specified and the seven failures are reported rather than hidden, and the theory states its conditioning (fixed partitions, f fixed) plainly. The one lapse is the Sc3 sentence calling a 0.012 difference "one and a half Monte Carlo standard errors" from an unpaired SE of 0.009; the paired SE is smaller and the difference is more significant than stated, so the text understates its own weakness.

**Clarity: 3/5 (weight 1).** The prose is precise but long: Sections 1 to 5 run to about 10,800 words, the compiled document to 73 pages, and the main text refers to "the earlier version" fourteen times, replicating a response-to-reviewers inside the paper. Paragraphs in Sections 3.2 and 3.5 carry 15 to 20 numbers each with no signposting. The single-arm paragraph quotes a ceiling of 160 at n_T = 200 while Table 4 prints prior ESS 172.3 (Gaussian) and 169.5 (survival) for the pooled model, above the ceiling, without explanation.

**Novelty: 4/5 (weight 2).** Every component is borrowed and the paper says so: the f + D g decomposition is BCF, the leaf prior is Hobbs' commensurate prior in spike-and-slab form, the ESS is Neuenschwander's ELIR. What is new is placing the spike-and-slab commensurate prior on the leaves of a discrepancy ensemble so that the partition decides where to borrow, the exact leaf marginals that make the sampler exact, the estimand-level ESS with a data-dependent ceiling, and the map P(|g(x)| < c | data) with its consistency result. The prior-cost argument for why one shared ensemble cannot do this is a useful piece of reasoning.

**Impact: 3/5 (weight 2).** The method is honest about its position: it beats BART-PP by 16% in RMSE when the sources agree and matches it under global conflict, but in the scenario it was built for it does not beat the free-offset model on the control surface (0.824 against 0.785 in R at delta_2 = 2) and loses at delta_2 = 1. The single-arm variant reduces to a commensurate AFT-BART with a flat map, and the application reaches the same point estimates as plain AFT-BART on the registry (0.97 and 0.90). The map and the ceiling are the deliverables a practitioner would use; the borrowing gain itself is modest.

**Performance: 3/5 (weight 2).** From Tables 1 to 3: Sc1 LRC-BART (100) RMSE 0.172 against BART-PP 0.205 and BART-NP 0.220, a clear win, with BART-CP and C-BART ahead at 0.142 and 0.150 because they pool. Sc2 RMSE 0.192 against 0.204 with bias 0.023 against 0.007. Sc3 bias 0.041 against 0.029 and RMSE 0.280 against 0.273, a loss; C-BART and BART-CP fail (0.23 and 0.98). Sc4 delta_2 = 1: RMSE 0.223 against 0.207, coverage 0.89, a loss on both X5 and X7; delta_2 = 2: ATE tied (0.202 against 0.205) but surface RMSE worse in both regions; delta_2 = 0.5: 0.200 against 0.203 at bias -0.043. Sc5: tied with BART-PP. Survival Sc1: median-ratio RMSE 0.312 against 0.345 and power 0.34 against 0.25, a win; Sc3: median-ratio bias 0.155 against 0.041 and RMSE 0.558 against 0.451, the largest loss in the paper, softened on the RMST ratio (0.146 against 0.137). Single-arm Sc2: every method biased by 0.4 or more and LRC-BART is not better than BART-CP. Against complete pooling and C-BART the method wins everywhere except Sc1; against BART-PP the record is two wins, three ties, three losses.

## Code audit (subagent, read-only)

| Item | Status | Evidence | Note |
|---|---|---|---|
| H_f=50, H_g=5, (0.5,3), nu_0=3, w~Beta(1,1) | PASS | 01-code/R/fit_lrcbart.R:38-43; lrcbart/R/lrcbart.R:116-120 | Harness passes no method_args, so wrapper defaults ran |
| Range-rule slab, k=2, (alpha_f,beta_f)=(0.95,2) | PASS | lrcbart.R:117-118,152-154 | |
| tau_0^2 truncated to (0,tau_1^2) | PASS | src/lrcbart.cpp:854 | |
| 1000 + 1000 iterations | PASS | comparators.R:66-67; fit_lrcbart.R:41 | |
| q=0.95, B=20, grid 141 on [1e-6,10], c_g by MC | PASS | lrcbart.R:358-367,407 | |
| sigma_1^2 from plain BART on RCT controls | PASS | lrcbart.R:385-393 | Stale comment at fit_lrcbart.R:67-69 says OLS/survreg |
| Map threshold c=0.5 | FLAG | fit_lrcbart.R:42; sections_1_2.tex:143 | "One third of sigma_1" is 0.467 for survival |
| DGP constants (mu, beta, f_0, alpha, delta, sigma, assignment, censoring) | PASS | dgp.R:26-28,46-55,69-98,142-149,222-223 | All match Section 3.1 and G.1 |
| Sc4/Sc5 shift sign | FLAG | dgp.R:204 adds +delta_2; G.1 says +delta_2, g=-delta_2; Section 3.1 says "shifted by -delta_2" | Section 3.1 wording contradicts code and G.1 |
| Single-arm n_T, sigma_1^2=sigma_2^2 | PASS | single_arm.R:31-33; lrcbart.R:399 | |
| Bias, SD, RMSE, coverage, power | PASS | metrics.R:49-59,121-128 | Match Section 3.1 |
| Surface RMSE profile set | FLAG | metrics.R:77-82 vs estimands.R:40 | Computed over all 300 RCT profiles, text says RCT control profiles |
| Seeds deterministic, replicate-paired | PASS | harness.R:40-51,94; dgp.R:168; run_*.R seed 2026 | Data seed 2026+r in every scenario |
| Timing claims | PASS | NOTES.md:45; results/timings.csv | |

No FAIL. The three FLAGs are wording, not computation, and cost half a point on Correctness.

## Top three strengths

1. The exactness of the sampler and the honesty of the theory. Every leaf marginal is closed-form, the verification reproduced a two-leaf joint posterior to Monte Carlo error, and the appendix says what is not proved (nothing about the posterior over partitions).
2. The ESS ceiling and the fractional ladder. Theorem 3(d) explains why an absolute target collapses under covariate shift, Sc2 shows it, and in the application the ceiling exposed a covariate-coding error the earlier analysis had missed.
3. The results are reported against the method's own pre-specified criteria, including the seven failures, and the Sc4 boundary leak is diagnosed (oracle f, every g prior, long chains) rather than tuned away.

## Top three weaknesses and fixes

1. The scenario the method was built for does not show it beating the free-offset model. Fix: add the pointwise ESS and the surface RMSE restricted to profiles more than 0.5 from the boundary, which the go/no-go report says separates the regions, and state the method's claim as "the map plus Sc1 efficiency" rather than local borrowing demonstrated on the surface.
2. No covariate-adaptive comparator. Fix: run CAHB (Jin et al. 2023) on Sc1, Sc4 (delta_2 = 1, 2) and Sc5, Gaussian only, R = 200; it is the closest parametric relative of the map and the referee will ask for it.
3. The survival Sc3 loss (median-ratio RMSE 0.558 against 0.451). Fix: report the RMST ratio as primary in Table 3 and the median ratio as secondary, matching the application, and give the paired SE of the difference so the reader can see that the RMST loss is within noise.

## Fixable issues (text and tables, no methodology)

- Section 3.1: "shifts the RWD surface by -delta_2" contradicts G.1 and dgp.R; write "+delta_2, so that g = -delta_2 outside R".
- Section 3.1: surface RMSE is over all N RCT profiles, not the RCT control profiles; correct the text or the code comment and Table 2 caption.
- Section 2.5: c = 0.5 is one third of sigma_1 for Gaussian only; say "0.5 on the outcome or log-time scale".
- Section 3.2, Sc4 delta_2 = 0.5: "g is -0.18 everywhere" against Table 2 values -0.15 and -0.18.
- Section 3.4: ceiling "160 at n_T = 200" against Table 4 prior ESS 172.3 and 169.5; explain why the pooled ESS_0 exceeds the ceiling or print the ceiling.
- Sc3 paragraph: replace the unpaired "one and a half Monte Carlo standard errors" with the paired SE.
- Build: appendix_new.tex redefines eq:model, eq:leafprior, eq:Mf, eq:Mg, eq:zpost, eq:ess0, eq:prrule; rename the appendix labels.
- Move the fourteen "earlier version" comparisons out of the main text into Web Appendix G and the cover letter; they are the revision history, not the paper.
- Cut Sections 1 to 5 by about a third toward the Biometrics length limit; the abstract at 248 words is at the limit.
- \author{} is empty; the comment block at the top of sections_1_2.tex about theorem numbering should not ship.
- Table 4 caption cites Web Table S6 while the appendix comment says tables run S1 to S5; confirm the numbering in the compiled PDF.
- Stale comment at fit_lrcbart.R:67-69 (OLS/survreg) should say plain BART.

## So what?

A trialist who borrows registry controls gets, for each patient profile, a posterior probability that the external data are compatible there and a local effective sample size, with a ceiling that says how much of the registry reaches the trial population at all. The efficiency gain over the safe free-offset model is 16% when the sources agree and nothing when they disagree on a region of moderate size, so the case for the method rests on the map and the calibration, not on the estimator.

## Verdict for Biometrics

Major revision. The methodology is sound, verified and clearly positioned, and the map with its ceiling is a contribution the journal's readers would use; but the paper is too long, still reads as a response letter in places, lacks the covariate-adaptive comparator its own motivation calls for, and needs its Sc4 claim reframed around what the tables show. With the fixable list applied and CAHB added, we would expect a B+ to A- on resubmission.
