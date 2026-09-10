# Reconciliation of the side-conversation full-paper review

The user supplied FULL-PAPER-REVIEW.md during final QA. It reviews the 82-page PDF created at 09:28:57 CDT, before this revision run. Its original major-revision verdict is retained as an assessment of that baseline. The review is a separate internal assistant assessment, not external peer review. The table below records the editing agent's reconciliation, not a new verdict from the side reviewer.

| Major finding | Current disposition | Candidate location and evidence |
|---|---|---|
| 1. Marginal ELIR | Addressed | Section 2.4 and Appendix D explicitly define a variance-based working precision summary, exact ELIR only for a normal prior/reference experiment. No predictive-consistency claim is retained. Primary Neuenschwander paper and deterministic mixture integration checked. |
| 2. Independent leaf indicators | Addressed | Section 2.4 gives exact weighted-leaf variance; D.1 averages all occupied leaf states. Tree-binomial shortcut restricted to one positive-weight leaf per tree. Realized/regional and application pointwise definitions traced to their producing functions. Mean partition coefficient remains an explicit calibration convention; Jensen difference demonstrated in checks. |
| 3. Calibration selection | Addressed | Appendix D states the rescaled block law, full-state chain assumptions, finite-grid selection, and no CLT/rate claim. Main text and implementation remark separate the ceiling cap and upper-grid fallback from the Pr-rule. Code unchanged. |
| 4. Ceiling | Addressed | Flat-prior one-tree anchor separated from finite regularization; no universal decrease under covariate shift. Truncated-law floor and finite-grid limits stated. Zero-weight empty leaves defined. |
| 5. G.6 | Addressed | Correct total configuration variance retains spike contributions. Optimizer derived for integer slab count; fixed-configuration density comparison distinguished from posterior probability. Clinical threshold/indicator inevitability removed from paper and slides. |
| 6. Support | Addressed | Section 2.2 and Appendix E distinguish population identification, shared-tree extrapolation, conditional empty-leaf law, and a globally single-arm likelihood. |
| 7. Single-arm survival | Addressed | Mean preservation restricted to linear control means. Section 3.4 and G.7 state nonlinear ratio means can move, with existing n_T=200 table examples. Section 4.3 explicitly gives PFS mean 0.97 to 1.00 under w=0.9. All numerical tables retained. |
| 8. Application | Addressed by limiting interpretation | Five baseline covariates and post-induction ASCT distinguished. Section 4.2 defines a descriptive model-based association conditional on recorded profiles, not a total causal effect. This review also prompted an explicit statement that equal residual variances govern counterfactual survival prediction, verified in lrcbart.cpp:451-463 and lrcbart.R:565-566. Changing ASCT adjustment or this assumption would require application reanalysis. No such reanalysis is claimed. |
| 9. CAHB | Addressed | Related work and G.0 correctly describe nonparametric kernel means/precision, current mean centered on historical mean. Transferred R=200 comparison caveats retained. Primary Jin paper checked. |

Writing: repeated introductory/tutorial prose shortened; main text now 6,993 words excluding the 227-word abstract by texcount text counts, before separate tables/captions. A/C/F development narration replaced with self-contained scientific scope. Double-log notation and stray comma corrected. Side review additionally prompted replacing residual claims of “unbiased” with empirical near-zero bias in the paper, appendix and slides. Author metadata remains unanswered and empty.

Presentation: hyperlink borders removed, preserving active links. Main numerical table and figure files are preserved byte-for-byte, including Table 2 and its continuation/regional panel; table selection, float restructuring and removal of the flat single-arm figure are optional presentation suggestions not applied in this scientific round. A diagnostic summary is supplied below for reviewers; the manuscript defines these quantities in Section 2.4, Section 2.5, and Appendix D without adding another main-text table.

| Diagnostic | Conditioning/averaging | Units and use |
|---|---|---|
| Compatibility B_c(x) | Posterior probability of absolute local discrepancy below prespecified c | Probability; local outcome agreement under the model |
| All-spike probability | Posterior event for the leaves reached by a profile | Probability of component labels; distinct from B_c |
| Calibration ESS_0 | All-spike reference; mean prior partition coefficient substituted before inversion; reciprocal variance averaged over untruncated scale law | Gaussian observation precision units; calibrates the spike-scale hyperparameter |
| Realized/regional ESS | Prior leaf variance evaluated at posterior partitions, indicators and scales; reciprocal variance averaged over draws | Working precision units; regional standardizations are not additive; not posterior variance reduction |
| Application pointwise ESS(x) | Stage-1 pointwise variance plus H_g times fitted mean spike variance, then inverted | Uncensored log-time reference precision; registry support under stated model, not RMST information |

These dispositions do not establish sampler convergence, validate the full implementation, resolve clinical design decisions, or certify submission readiness. The stored code-impact note records what would require new calibration or results.
