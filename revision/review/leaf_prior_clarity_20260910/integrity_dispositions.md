# Integrity disposition

The user explicitly authorized the variance index to depend on z after the initial prose clarification. The later authorization extends equation (2) and matching formulas; the initial instruction to preserve equation (2) is therefore superseded within this bounded notation change.

Checks use `check_latex_rewrite.sh --allow-structure-change --allow-authorized-content-change` with prose/slides modes. Main and appendix checks return 1 solely because their numbered-equation text changed. Slide checks return 0; its unnumbered mathematics still requires the recorded manual diff inspection.

- Main equation (2): binary-selector mixture replaced by the equivalent conditional normal with variance tau_{1-z_hl}^2, and scales/w made explicit in conditioning.
- Main unnumbered displays: added the definition of g(x), all-slab f_1 conditional, and general f_1 conditional; the preexisting all-spike display remains. Conditional variance is the sum of independent reached-leaf variances. Replacing one term tau_0^2 by tau_1^2 gives the stated variance increment.
- Main conditional theta draw: custom tau_z notation replaced by tau_{1-z}.
- Appendix equations eqA:leafprior, eqA:Mg, eq:leafcond and eq:linform: same index mapping, with a two-line layout for eq:leafcond. The rest of each formula is unchanged. Inline proof occurrences use the same mapping, and the joint density is displayed for line fit. Continuous theta is integrated out, replacing the imprecise phrase summing out.
- Slide 8: same indexed conditional normal, explicit theta definition and index-to-component mapping.
- Number-count deltas arise from index 1-z, component labels and explanatory mathematical displays; no empirical number, hyperparameter value or reported result changed. Structure deltas are the added unnumbered displays, aligned environment and samepage wrapper.

Citations, labels, sections, references, graphics and protected table/figure blocks pass. Main wrapper, Sections 3-5 and Appendix G match the local follow-up baseline byte for byte. All three source diffs have been inspected. Hyperpriors, truncation, scales, existing scientific claims and result files are preserved.

Curated source diffs use zero context to keep the archive whitespace-clean; full-context diffs and original raw build logs are preserved in the local follow-up run. Curated compiler/audit text has trailing whitespace normalized.
