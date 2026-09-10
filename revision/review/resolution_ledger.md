# Technical resolution ledger

- H1: Approved moves completed; source blocks preserved in baseline and cut_log.
- H2: CAHB primary rows generated from R=200 CSV; all existing rows retain R=500. CAHB map/ESS columns suppressed in Table 2 because definitions differ; separate appendix panel gives precision shares.
- H3: Regional ESS panel uses R=200 and all-RCT regional standardization. Control-profile and far-region panel generated from paired.csv. All-profile definition retained.
- H4: Abstract, simulation and Discussion distinguish posterior compatibility, realized regional ESS and surface performance. No general ranking is claimed for the transferred CAHB analysis.
- C1: Main theorem's truncated-prior zero-limit/convexity claim contradicted Appendix D. Main theorem now restricts those statements to the untruncated calibration law and states the positive truncated floor. Protected appendix unchanged apart from authorized implementation remark.
- C2: Main theorem block asymptotics were missing the appendix's fixed-length/increasing-length distinction and rescaling conditions. Summary qualified to the appendix's regimes.
- C3: Student-t leaf marginal qualified to the untruncated scale law; sampler truncation explicitly distinguished.
- C4: w^H_g describes leaves covering a profile, not every leaf in an ensemble; wording corrected.
- C5: Application ESS calibration concerns the control log-time mean, not RMST; corrected to match calculation.
- C6: Regional shift 1 at residual SD 1.5 is two thirds of a residual SD, not one; corrected. Survival ratio biases no longer called an unbiased ATE.
- C7: Unsupported claims of universal BART-PP safety, no agreement gains, a proved sample-size limit and additive multi-source ceilings removed or narrowed. Failed tuning configurations remain documented.
- C8: Historical application comparisons 1.24/1.25 versus 1.21/1.18 are not equal to rounding. Wording corrected; observed tail reduction is not called a bounded tail.

Remaining scientific review items (not resolved by prose): the prior-averaged binomial ESS formula versus independent per-leaf indicators on nontrivial partitions; exact ELIR interpretation for a nonnormal marginal prior; source-level identification outside trial support; formal validity of the shared-leaf prior-cost threshold in G.6. Preserve these mathematical objects pending scientific review. The post_technical directory is a prose-validation baseline, not a certification that these issues are closed.

## Checker exception mapping

The original-to-post_technical comparisons are intentionally not ordinary prose checks. Their full output is preserved.

| File | Hard-check categories | Mapped changes |
|---|---|---|
| Sections 1-2 | numbers, citations, references; inline-math audit | H1 moves remove commentary/calibration tokens and move rockova2020bart to G.3; C1-C4 correct the theorem and distributional/profile qualifiers. No displayed equation or label was removed. |
| Sections 3-5 | numbers, citations, references; inline-math audit | H1 moves coding details to G.5; H2 adds Jin2023CAHB and its R=200 results; H3-H4 add regional ESS and paired surface summaries; C5-C8 correct the estimand, shift units and unsupported interpretations. |
| Appendix G | numbers, citations, references; inline-math audit | H1 receives moved commentary and coding details; H2 receives CAHB transfer and variants; H3 receives control-profile and far-region summaries. |

The final prose comparisons use post_technical and pass without allowance flags. Number-word deltas are rhetorical counts or ordinals (one/first/second/three) removed when restructuring sentences and lists; scientific counts and thresholds are unchanged. The Sc3 explanation was explicitly qualified as a possible explanation because the conditional leaf theorem does not establish a decomposition of full-ensemble bias (C7). G.6 remains review-retain.

Appendix A-F preservation check: removing the newly inserted calibration implementation remark yields a byte-identical copy of the original appendix_new.tex. Wrapper-only typography changes use microtype and an emergency-stretch setting; scientific text is unchanged.

Layout exceptions: long main tables continue across pages with repeated headings and intact scenario blocks; supplementary hyperlink anchors now have a separate prefix. These changes do not alter table cells, labels or citation targets.
