# Summary

Round-2 applied/contextual resolution review of the candidate manuscript assembled by `candidate/05-writing/main_draft_critical_final.tex`, its main text and Appendices A–G, and the candidate slide deck. I consulted my own round-1 review, not the other reviewer or combined report, and made no manuscript, code, table, or slide edits. This is an independent AI review, not external human approval.

All five prior major comments are resolved. The candidate now distinguishes a working variance-based information summary from exact marginal ELIR, treats G.6 as a constrained prior-density illustration, accurately positions CAHB as kernel based and nonparametric, distinguishes off-support prediction from identification, and presents the application descriptively with post-induction ASCT adjustment. It also corrects the nonlinear survival-ratio mean claims using existing results. Three short residual corrections remain; none requires a new analysis or alteration of results.

# What the paper does well

- The external-centered model explanation and diagram are preserved.
- The revised calibration exposition specifies squared leaf proportions, indicator aggregation, partition averaging before inversion, the untruncated calibration convention, fitted truncated prior, and implementation overrides. The application map's plug-in scale and the supplementary prior precision average are distinguished.
- The deck now gives the regional surface loss away from the boundary, the distinct regional ESS, correct residual-SD units, and CAHB's paired `R=200` transfer qualification.
- The application labels ASCT after induction correctly and avoids an efficacy conclusion; the earlier-tail comparison no longer assigns a cause that the unavailable sampler cannot verify.
- The main text is materially shorter without sacrificing the central support, calibration, comparator, or empirical-loss caveats.

# Resolution check

- **M1 — Resolved.** `sections_1_2.tex:44,93` motivates separate partition/shrinkage priors and explicitly says the calculation “does not establish a detection threshold ... or rule out regional borrowing with shared leaves.” G.6, `web_appendix_G.tex:127–141`, retains integer slab counts and all remaining spike contributions, gives the constrained minimum, and distinguishes joint-density modes from posterior probabilities. The deck's “Why a separate discrepancy ensemble,” lines 216–226, follows this limited interpretation.
- **M2 — Resolved.** `sections_1_2.tex:49,51` calls CAHB the closest kernel-local relative, expands the name correctly, and recognizes its nonparametric kernel estimates. G.2, lines 24–30, retains all transfer choices, distinct precision-share definition, and variant results. The deck identifies the transfer and first 200 paired replicates at lines 319 and 368. One residual current/historical wording error in G.0 is listed below.
- **M3 — Resolved.** `sections_1_2.tex:93` distinguishes source-support identification, extrapolation, empty-leaf conditional priors, and information transmitted through shared leaves/hyperparameters. Appendix E, lines 477–479, explicitly says off-support nonidentification “does not imply that its fitted posterior law ... equals the prior,” while preserving the true no-controls-anywhere result.
- **M4 — Resolved.** `sections_3_5.tex:113` distinguishes five baseline covariates from ASCT after induction; line 122 calls the comparison descriptive and states that it does not identify a total causal effect. Results at line 129 are conditional descriptive comparisons. G.5's cohort caption and the application slides are synchronized; no causal reanalysis is claimed.
- **M5 — Resolved.** `sections_3_5.tex:101,131` distinguishes Gaussian centering from nonlinear survival-ratio means and reports the observed PFS sensitivity change from 0.97 to 1.00. G.7, lines 151–153, gives the existing median-ratio bias contrasts 0.076/0.108 and 1.492/0.108 rather than claiming invariance. The single-arm slide, line 411, makes the same distinction.

# Major comments

None.

# Minor comments

## m1. One G.0 sentence still puts the CAHB prior on the wrong mean

1. **Location:** `candidate/05-writing/web_appendix_G.tex:15`: “places a commensurate prior on the historical control mean function.”
2. **Issue and consequence:** CAHB's prior is on the current control mean, centered on the historical mean function. G.2 line 24 already states this correctly; G.0 reverses their roles. The primary-paper verification from round 1 remains sufficient: [Jin et al. (2023), Section 2](https://jinhuaqing.github.io/files/2023_SIM_CAHB.pdf).
3. **Concrete fix:** Replace this clause by “places a commensurate prior on the current control mean function, centered on the historical mean function.” No other literature or comparator change is needed.

## m2. Retain the fixed-scale, large-count qualification when interpreting abstention

1. **Location:** `candidate/05-writing/web_appendix_G.tex:57`: “a shift is detected only beyond ...” and “abstention under a shift of several spike scales is near certain”; related numerical extrapolation in `candidate/03-theory/appendix_new.tex:313`.
2. **Issue and consequence:** Theorem A2 gives a conditional detection probability and a fixed-scale large-count limit. A small calibrated prior scale alone does not establish near-certain finite-sample abstention for the learned-scale, learned-partition ensemble. A probability threshold is also not a strict impossibility boundary for smaller shifts. These two sentences read more strongly than the now-correct main text and surrounding caveats.
3. **Concrete fix:** Say the conditional half-detection threshold depends on the spike scale and leaf standard error, and that the limiting spike probability can be small for discrepancies large relative to a fixed spike scale. Replace the learned-prior inference with one sentence that finite-sample ensemble behavior is evaluated empirically. Keep all theorem formulas, numerical illustrations, and simulation results.

## m3. Fix two small literal remnants

1. **Location:** `candidate/05-writing/sections_1_2.tex:103`: “A g leaf uses residuals ... use only RCT controls”; `candidate/06-slides/main_lrcbart.tex:410`: a capped fit “reduces to pooling.”
2. **Issue and consequence:** The first is a duplicated verb introduced during shortening. The second describes the finite smallest-grid spike scale as an exact zero-scale limit, although the single-arm slide assumes `w=1` and Appendix D correctly calls the result an approximation.
3. **Concrete fix:** Use “A g leaf uses residuals ... from RCT controls only.” Change the slide to “approximates pooling with w=1.” Neither correction changes results or the slide structure.

# Suggestions to enhance the paper

No further structural or literature agenda is needed this round. Finish the bounded corrections above and the planned build, visual, redline, and numerical-lineage checks. Any future move to exact marginal ELIR, altered calibration law, or a causal application estimand would require separately specified analyses; the current prose appropriately avoids claiming those changes.

# Source-dependent flags

- The author field remains empty. Final author order and affiliations require author-supplied metadata; no inference from deck credits is warranted.
- The clinical interpretation is now appropriately descriptive. Independent clinical qualification of endpoint/covariate definitions and the external-control comparison remains outside this AI manuscript review.
- The former MAP-AFT-BART sampler remains unavailable. This limitation is now disclosed and no longer supports an asserted causal attribution of the earlier tails; the source flag is handled by narrowing the claim, not by reproducing that analysis.
- This applied review assesses the coherence of the revised working-ESS interpretation. It is not a new empirical calibration validation, sampler convergence audit, or external mathematical approval.

# Verdict

The prior blocking applied/contextual issues have been resolved without replacing results. The remaining corrections are local and do not require new analyses. The revised manuscript can proceed through artifact QA with these minor edits and the stated author/source boundaries preserved.

VERDICT: minor revision
MAJOR COUNT: 0

What bothered me most: a few isolated appendix sentences still retain the stronger interpretation that the revised main text has successfully removed.
