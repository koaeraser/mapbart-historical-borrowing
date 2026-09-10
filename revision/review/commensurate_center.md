# Commensurate-prior center clarification

User-requested presentation correction, September 10, 2026.

The paper's Section 2.2 and slides 5–8 now identify f(x) as the external-control response mean and show the commensurate distribution of the RCT-control mean centered at f(x). Conditional on the partitions, scale and all reached leaves being in the spike, its variance is H_g tau_0^2. The discrepancy is then introduced as f_1-f, explaining its zero-centered prior. The paper also gives the general conditional variance as the sum of reached-leaf variances and distinguishes latent-mean uncertainty from residual variance.

The numbered leaf-prior equation (2), sampler, numerical results, tables and protected theory are unchanged. One explanatory slide was added; the tree diagram and model/leaf-prior wording were updated. The current paper is 82 pages, the cumulative tracked paper 90 pages, and the deck 25 slides.

Validation: paper and tracked-paper builds have no warnings or overfull boxes. Slides have no overfull boxes or unresolved references; the existing Metropolis/pdfLaTeX and font-substitution warnings remain. Changed pages/slides were inspected at full size and the complete clean paper/deck at overview scale. The tracked model pages were inspected separately. A local sloppypar in the tracked source accommodates deletion/addition markup without altering the clean source.

The prose checker passes citations, labels, sections, numbered equations, protected blocks, graphics and structure. Its additional reference to equation (1) is an audited exception: the new sentence distinguishes response-mean variance from observation residual variance. Added formula/inline-math tokens express the authorized centering explanation. Slide integrity checks pass citations, labels, sections, numbered equations, protected blocks, graphics and references; the new frame and explanatory formulas are authorized additions.

Remaining scientific-review questions recorded in the existing ledger are unaffected.
