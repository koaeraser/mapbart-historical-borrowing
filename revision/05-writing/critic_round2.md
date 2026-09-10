# Round 2: internal consistency and editorial review

Verdict: TEXT FIXES APPLIED; SCIENTIFIC REVIEW STILL REQUIRED.

This is a review by the same Codex agent that edited the manuscript, not an independent paper-critic assessment. It must not be represented as independent scientific approval.

## Findings resolved in the text

1. **Major: calibration theorem disagreed with its proof.** The original main Theorem 3(a) asserted a zero large-scale ESS limit and convexity under truncation. Appendix D explicitly gives a positive floor and claims convexity only without truncation. The main statement now makes this distinction. Its Monte Carlo summary also refers to the appendix's separate block-number and block-length regimes and rescaling conditions.
2. **Major: map and ESS were conflated.** The abstract and Methods described a posterior compatibility probability as an amount of information used. The text now distinguishes compatibility, pointwise prior information and realized regional ESS. Table 2 adds the latter at R=200, with nonadditivity and regional standardization stated.
3. **Major: CAHB diagnostics are not LRC-BART diagnostics.** CAHB's map is 1-1/R(x), a precision share; its effective size is CECSS excess. These values are not printed under LRC-BART's B_c or realized-ESS headings. Appendix G.2 states the five transfer decisions, diagnostic-seed overlap, interval construction, and p10/literal variants. Existing comparisons use matched R=200 values in the prose.
4. **Major: Sc4 interpretation exceeded the evidence.** The regional surface loss remains after boundary exclusion. The text reports its paired SE and distinguishes all-RCT from control-profile standardization. Failed tuning is no longer described as proof of an unavoidable sample-size limit.
5. **Minor: distributional and estimand wording.** The Student-t leaf marginal is qualified to the untruncated scale prior. The w^H_g probability concerns leaves covering a profile, not every leaf of every tree. Application calibration concerns the standardized control log-time mean, not RMST. A shift of 1 at residual SD 1.5 is two thirds of a residual SD. Survival ratios are no longer called an unbiased ATE.
6. **Minor: application comparisons.** Values 1.24/1.25 versus 1.21/1.18 do not agree to rounding. That description was removed; finite observed tails are no longer called bounded tails.

## Unresolved scientific questions

- The prior-averaged ESS formula uses a Binomial(H_g,w) sum. The model assigns independent indicators per leaf, and an estimand can span several leaves of each tree. The formula's exactness for nontrivial partitions needs a derivation against the implementation's leaf-weighted discrepancy variance.
- The identification of inverse prior variance with exact ELIR needs qualification for the BART/mixed marginal prior. Conditional normal information and information in a marginal mixture are different quantities.
- The statement that g is a prior draw outside trial support needs to distinguish absence of a local likelihood contribution from posterior information transmitted through shared tree leaves and hyperparameters.
- The G.6 prior-cost calculation needs a separate scientific check of its optimization constraints and its transition from prior costs to indicator behavior. It is marked review-retain and was not rewritten.

These are not resolved by a clean build or a token-invariance check. The protected appendix was retained except for the authorized Appendix D implementation remark. The manuscript remains a working candidate; the author field is still empty.

## Evidence

See humanize_logs/codex_20260910/resolution_ledger.md, numerical_audit.md, the strict *.integrity.log files, and the immutable pre-edit and post-technical baselines. Thirty new numerical claims were checked against scenario/method/metric-specific CSV cells; all original numerical table rows were retained. No simulation or application results were modified.
