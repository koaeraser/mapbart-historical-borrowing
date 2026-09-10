# Next revision: scope and authorization

The user requested the next revision round under HANDOFF-NEXT-REVISION.md. The current clean build is up to date. Original and restart candidates remain immutable; all editing starts in candidate/. Independent reviewers do not edit manuscript sources. The critical-revisor review/revision cycle is used with verification following integration. The requested deliverables retain the handoff paths after verified promotion.

| Finding | Allowed objects | Evidence for closure | Fixed invariants |
|---|---|---|---|
| R1 ESS aggregation | Section 2.4, Appendix D, slide ESS formulas | Leaf-weighted variance against lrcbart.R; analytic checks | Existing simulation and application values, sampler |
| R2 ELIR | ESS definitions/interpretation and connected claims | Primary ELIR paper; conditional versus marginal information derivation; calibration source | Existing calibration results identified under implemented rule |
| R3 Support | Section 2.2, Appendix E if necessary, slides | Model likelihood factorization and tree/hyperparameter sharing | External-mean centering |
| R4 Prior cost | G.6 and connected model/novelty claims/slides | Constrained optimization and distinction from posterior selection | No unrun empirical claims |
| R5 Reconciliation | Full manuscript/deck caveats, unsupported guarantees | Reviewed baseline and producing code/CSV/source | Tables, plots, sample counts, comparator definitions |
| R6 Length and authors | Main-text prose, verified author metadata | texcount; user answer for author order | Scientific caveats and distinctions |

Any additional supported correctness finding from the independent reviews will be explicitly mapped before editing. Appendix A-F changes require a scientific reason in the ledger. No full simulations are authorized by implication from wording corrections. Existing upload authorization is retained as stated in the handoff; no separate message to Yunxuan will be sent.

## Accepted independent-review findings before integration

- Reviewer 1 M1-M5 and m1-m4: working ESS, leaf aggregation, block-law/cap/anchor correction, population support, constrained prior cost, pointwise convention, interior w, nonlinear means, invariant-kernel wording. Direct proofs/code and deterministic checks support these changes.
- Reviewer 2 M1-M5 and m1-m4,m6: same scope issues plus verified nonparametric CAHB classification/name, ASCT after induction and descriptive interpretation, actual cohort masks, coding extrapolation, log-time notation and log-RMST diagnostics. Producing source/code/CSV support these changes.
- Both reviewers' prose/history comments: targeted main-text shortening and removal of obsolete proof references are allowed, with Appendix A-F changes confined to connected scientific corrections. Author metadata remains pending user input.
- Original tails may be reported, but their cause cannot be isolated with the earlier sampler unavailable. Narrowing that attribution is required; no old result is altered.
