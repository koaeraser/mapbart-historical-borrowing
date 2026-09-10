# Summary

Round-2 independent AI methodologist resolution review of the candidate assembled by `candidate/05-writing/main_draft_critical_final.tex`. I used only my own round-1 review, the candidate manuscript and its source changes, and the algebra/code evidence in `checks/mathematical_checks.txt`; I did not read the other reviewer or combined review. The five previous major objections are resolved in substance. The remaining points are small formulation and editorial corrections. This review does not certify full-model convergence, operating characteristics, or clinical/causal qualification.

# What the paper does well

The revision now defines the implemented information summary, gives the actual leaf-proportion formula used in realized ESS, and states the calibration approximations and overrides. The corrected Appendix D proof is self-contained on its mathematical core. G.6 now compares the intended constrained prior objective and explicitly avoids posterior claims. The application correctly distinguishes a descriptive adjusted comparison from a total causal treatment effect, and nonlinear survival means are no longer claimed to remain unchanged under prior noise.

# Resolution check

- **M1 — Resolved.** `candidate/05-writing/sections_1_2.tex:122-139` explicitly calls the quantities working precision summaries, states the normal reference experiment, and says they are not exact marginal ELIR. `candidate/03-theory/appendix_new.tex:382-409` distinguishes marginal log-density curvature, component averaging, the mean-partition substitution and inverse marginal variance. Uncensored log-time reference units are stated. The old predictive-consistency implication is removed.
- **M2 — Resolved.** Main text lines 124-137 and Appendix D lines 390-409 define `V_g=sum a_hl^2{z_hl t+(1-z_hl)tau1^2}` and average independent states over L_+ leaves. The pointwise/stump binomial reduction is correctly restricted. Main line 149 and Appendix lines 464-466 describe the exact squared-leaf-proportion realized output as a prior-variance diagnostic evaluated at posterior states. Existing numerical outputs can be retained. The unused tree-state `lrc_ess_prior` helper remains a documented implementation limitation, not the source of the reported results.
- **M3 — Resolved in substance.** Appendix D lines 429-454 define the correct rescaled law `gamma U`, impose threshold continuity and finite-grid separation/nonempty selection conditions, and use a valid ergodic sandwich argument. The unsupported CLT and rate are removed. Lines 461-466 distinguish finite inner Monte Carlo from the exact-integral theorem and identify the cap and upper-grid fallback as separate branches. The flat-prior single-tree restriction, finite-prior information and zero-RWD-leaf issue are now explicit. The range-rule identity correctly shows that increasing H_f alone does not remove prior information. Minor m1 below refines the Markov-chain assumption to refer to the full sampler state.
- **M4 — Resolved.** Main line 93 and Appendix E lines 472-487 separate unrestricted population identification from posterior extrapolation through partitions. The text no longer claims that every unsupported profile is a prior draw; the leaf-wise conditional-prior result and single-arm special case are correctly separated. The single-arm variance is conditional on shared scales and w. Sc2 is described using limited overlap rather than a population support hole.
- **M5 — Resolved.** G.6 lines 127-141 minimizes the joint prior cost over all contributions with `S_m=(H-m)tau^2+m lambda^2`, restricts m to integers, and gives the correct comparison with the all-spike mode. It removes the unsupported numerical threshold and indicator-flipping claims. The pooled Gaussian likelihood loss now uses the variance of the between-source mean difference. Its scope as a conditional prior illustration is clear.

Prior minor comments m1-m4 are also resolved in substance: the application plug-in and supplementary integrated pointwise definitions are distinguished; main conditional theorems require interior w; zero-mean preservation is restricted to linear means with nonlinear examples supplied; Appendix A now says the composite kernel preserves the joint conditional. Prior minor m5 is only partly resolved because internal production references remain outside Appendix D; see m4 below.

# Major comments

None.

# Minor comments

## m1. Put the Markov assumption on the full Stage-1 state

**Location.** `candidate/03-theory/appendix_new.tex:429`: “Let the standardized Stage-1 draws form a ... Harris-ergodic chain.”

**Issue and consequence.** A scalar function of the full BART state need not itself be Markov, even though the full sampler is a Markov chain. The proof only needs the standardized surface to be a square-integrable measurable function of the stationary ergodic state and its nonoverlapping blocks.

**Concrete fix.** Write: “Let the Stage-1 state chain be stationary, aperiodic and Harris ergodic, and let its standardized surface draws be a measurable function with finite second moment and variance V_*>0.” Keep the rest of the proof; the blocking and sandwich argument then apply directly to the relevant sequence. This is an assumption statement, not a claim that the implemented BART chain has been proved to satisfy it.

## m2. Specify the zero-weight convention in the anchor sum

**Location.** Appendix D lines 436-440 and main Theorem 3(d).

**Issue and consequence.** The allowed case a_l=pi_l=0 produces a formal 0/0 in the displayed sum despite contributing no information. The intended condition excludes pi_l=0 only when a_l>0.

**Concrete fix.** Define `sum a_l^2/pi_l` over leaves with pi_l>0, assigning zero contribution to a_l=pi_l=0. Alternatively state this convention once. No numerical or theoretical change is needed.

## m3. Repair a duplicated verb in the posterior-computation paragraph

**Location.** Main text line 103: “A g leaf uses residuals ... use only RCT controls in that leaf”.

**Issue and consequence.** The shortening introduced a grammatical break in the description of the likelihood inputs.

**Concrete fix.** Use “A g leaf uses only RCT controls in that leaf, with partial residuals ..., count ..., and mean ...”. Preserve the formulas.

## m4. Finish removing internal production references from the supplement

**Location.** Appendix C lines 312-316 and 378 (“the specification”, “verification report, check 6c”); Appendix F lines 490-492 (superseded theorem mapping); analogous old-version narration in Appendix A's remark.

**Issue and consequence.** These unresolved parts of prior m5 make the supplement partly a revision history and depend on a report that a scientific reader does not have. The current Appendix D no longer has this problem.

**Concrete fix.** Replace the report reference with the actual current simulation section, state the fixed-scale/nonvanishing-indicator limitation directly, and put old-to-new theorem accounting in the ledger. Appendix F can retain its section slot with a short self-contained explanation of why conditioning on the global scale permits exact Gaussian-mixture leaf integration. Preserve the substantive conditional guarantees and gaps.

# Suggestions to enhance the paper

No new scientific agenda is needed in this round. Complete the four small corrections and retain the stated separation between mathematical results, working information diagnostics and empirical performance.

# Source-dependent flags

- The original MAP-AFT-BART sampler remains unavailable. G.5 now appropriately reports observed tail differences without claiming an isolated causal attribution.
- Code/results were preserved; no operating-characteristic reruns occurred. The deterministic checks verify bounded algebraic and producing-code identities, not complete sampler correctness or mixing.
- The unused exported tree-wide-state prior-ESS helper should remain clearly documented in the implementation-impact record. A future replacement of calibration with exact partition averaging or exact marginal ELIR requires a new calibrated empirical analysis.
- Author metadata and full clinical/causal qualification are outside this methodologist resolution review.

# Verdict

The revised manuscript resolves the material mathematical and interpretive objections from this review. The remaining edits are minor and do not require changing fits or numerical outputs. Acceptance here concerns the bounded manuscript revision; it is not independent external approval or a submission-readiness certification.

VERDICT: minor revision
MAJOR COUNT: 0

What bothered me most: a few old production references survive after the scientific core has become substantially more self-contained.
