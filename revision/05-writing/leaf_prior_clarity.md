# Leaf-prior clarification and indicator-dependent variance

September 10, 2026. This follow-up implements the user's requests to define theta_hl, explain z=0 for f_1, and make the variance index depend explicitly on the indicator. It follows the completed scientific revision; the earlier review reports remain assessments of their archived versions.

Equation (2) now states `theta_hl | z_hl, tau_0^2, tau_1^2 ~ N(0, tau_{1-z_hl}^2)`, with `z_hl | w ~ Bernoulli(w)`. The index is `1-z` because the existing convention is z=1 for the spike tau_0^2 and z=0 for the slab tau_1^2. This is algebraically identical to the previous binary-selector form. The notation is reconciled in the main conditional leaf draw, Appendix A prior, marginalization and conditional formulas, Appendix B covariance matrix, and slide 8.

The model description defines theta_hl as the discrepancy contribution of leaf l in tree h and l_h(x) as the reached leaf. It displays g(x) as the sum of reached-leaf contributions. Conditional on f, partitions, indicators and scales, f_1(x) has mean f(x) and variance equal to the sum of tau_{1-z}^2 over the reached leaves. Switching one leaf from spike to slab increases this variance by tau_1^2-tau_0^2. All-spike and all-slab configurations have variances H_g tau_0^2 and H_g tau_1^2, respectively. Both remain centered at f(x).

Validation:

- Clean paper: 79 pages; cumulative tracked paper: 101 pages; synchronized deck: 25 slides. The main text has 7,188 texcount text words including the 227-word abstract, or 6,961 excluding it. This uses the same two main section files as the prior count and excludes headings, captions, mathematical-expression counts and separate tables.
- Clean and cumulative tracked builds succeed with no LaTeX warnings, unresolved references/citations or overfull boxes. The slide build retains existing Metropolis/font substitution warnings, with no overfull boxes or unresolved references/citations.
- Targeted rendered-page inspection covers the model displays, Appendix A leaf marginal/conditional/proof, Appendix B covariance, and slide 8. The longer index required a two-line appendix conditional and a displayed proof density. Both clean and tracked versions fit their margins.
- Citations, labels, section structure, references, graphics and protected table/figure blocks pass the rewrite checks. The equation-preservation checks correctly flag the authorized indexed formulations; the accompanying integrity disposition maps every mathematical change. There is no blanket prose-only pass.
- The main wrapper, Sections 3-5 and Appendix G are byte-identical to the follow-up baseline. Code, numerical tables, figures, results and portable distribution paths were not edited. No simulations, application fits or additional independent scientific reviews were run for this notation follow-up.

Local baseline and detailed evidence: `05-writing/humanize_logs/leaf_prior_clarity_20260910_103907/`. The cumulative redline still compares the original pre-Codex baseline, using the existing table-aware generator from `05-writing/critical_revisor_logs/run_20260910_094551/build_tracked.py`. Prior scientific review snapshots are preserved. The distribution baseline is commit `594877e1e5018276e8f2f69de148a3eb72327768`; the curated follow-up record is `revision/review/leaf_prior_clarity_20260910/` in the existing PR.
