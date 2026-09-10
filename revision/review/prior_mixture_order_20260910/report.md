# Prior mixtures in model order

September 10, 2026. Section 2.2 now presents two prior displays in sequence: the two-component leaf mixture in equation (2), followed by the induced pointwise mixture for f_1(x). The equivalent indicator representation with variance tau_{1-z_hl}^2 remains in prose. The all-spike and all-slab components are explained below the induced mixture, with their variances and common center f(x).

The induced prior has H_g+1 components indexed by the number m of reached spike leaves. Conditional on f, partitions, w and the scales, its weights are choose(H_g,m) w^m (1-w)^(H_g-m), and its component variances are m tau_0^2+(H_g-m) tau_1^2. This follows from the independent reached-leaf indicators at a fixed x; it does not change the standardized-mean weighted-leaf calculation.

Only the prior exposition in Section 2.2 changed. Hyperpriors, scales, truncation, the Student-t qualification, every other manuscript section, appendix, slides, code, tables, figures and results are preserved. Author names and affiliations remain intentionally blank. No simulation or application fits were rerun.

The clean paper remains 79 pages and the cumulative tracked paper 101 pages. Main-text texcount is 7,159 words including the 227-word abstract, or 6,932 excluding it, using the same two main section files and exclusions as the preceding count. Both builds have no warnings, unresolved references/citations or overfull boxes. Clean pages 4-6 and tracked pages 7-9 were inspected at full-page scale; complete clean/tracked overviews were inspected. The two prior displays fit on clean page 5. The cumulative redline retains the original pre-Codex baseline and existing table-aware generator.

The rewrite checker passes citations, labels, sections, references, graphics and protected blocks. Its sole hard failure is the authorized equation (2) change. The integrity disposition records the mixture derivation and every mathematical delta. The source outside the bounded prior passage is byte-identical to the follow-up baseline.

Local baseline and evidence: `05-writing/humanize_logs/prior_mixture_order_20260910_123054/`. Distribution baseline: `e33430e85044581f2fdb125069727f12f0511bb1`. Earlier review reports, baselines and manifests remain historical records. Current distribution evidence is in `revision/review/prior_mixture_order_20260910/`.
