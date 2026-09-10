# Round 2 integration

Both reviewers report that their five original Major comments are resolved. Remaining pinpoint comments were integrated as follows:

| Finding | Disposition | Change and evidence |
|---|---|---|
| Method reviewer: projected chain | Applied | Appendix D(c) assumes a full-state Harris-ergodic chain and a square-integrable standardized-surface function. Main theorem summary reconciled. A measurable projection need not itself be Markov. |
| Method reviewer: empty zero-weight leaves | Applied | Defined zero contribution when both target and RWD leaf proportions vanish, in the main theorem and Appendix D. |
| Both: duplicate verb | Applied | Corrected the g-leaf residual sentence. |
| Context reviewer: CAHB prior target | Applied | G.0 prior is on current mean centered at historical mean, matching the primary paper and G.2. |
| Context reviewer: G.3 abstention | Applied | Qualified the large-count fixed-scale limit in G.3 and Appendix C; no learned-scale finite-sample guarantee. Replaced the strict detection statement with a conditional half-detection threshold. |
| Context reviewer: pooling wording | Applied | Single-arm slide says approximates pooling with w=1 at a finite minimum scale. |
| Method reviewer: obsolete proof references | Applied | Removed specification/check-report dependencies in A/C and replaced F with the scope of posterior invariance, separating invariance from independent draws/mixing. Baseline preserves the historical account. No theorem B/C equations or current numerical results changed; obsolete rough threshold .7 was removed from historical parenthesis while actual d_half=.54 remains. |

The A/C/F edits are documented scientific scope and self-containment corrections under the handoff exception, not a broad rewrite of protected proofs. Author metadata remains pending. No full-model operating-characteristic claim is newly certified.

## Side-conversation review received during QA

See `../side_review_reconciliation.md` for all nine major findings. The reviewed version was the pre-run 82-page baseline. Residual empirical-unbiasedness wording and the single-arm residual-variance assumption were corrected after checking the producing code. Active hyperlink borders were hidden as a low-risk presentation fix. This is additional integration of the supplied review, not a third commissioned review round.
