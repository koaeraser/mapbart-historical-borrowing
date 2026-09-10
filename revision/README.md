# LRC-BART manuscript and code revision

This directory contains the September 10, 2026 manuscript candidate, its implementation, and the summary results used to produce the tables. The original repository files remain available alongside it.

Start with the [paper](05-writing/main_draft.pdf) and [review ledger](review/resolution_ledger.md). The paper includes the CAHB comparison, regional effective sample sizes, and the Sc4 boundary analysis, with a revised exposition throughout. This is a working candidate for scientific review; the editing reports are not independent scientific approval.

## Review requested

Yunxuan, please review the manuscript and code, particularly:

- The prior-averaged ESS calculation when an estimand spans several independently selected leaves.
- The interpretation of ELIR and inverse-variance information for the BART mixture prior.
- Identification outside trial support and information shared through leaves and hyperparameters.
- The shared-leaf prior-cost argument in Appendix G.6, which was retained pending scientific review.
- The CAHB comparison, application definitions, and agreement between the implementation and the reported results.

The author field still needs completion. The main text is approximately 7,500 words, excluding the abstract, and needs further shortening if the 7,000-word target applies.

## Contents

- `01-code/`: R/C++ sampler package, comparators, scenarios, and tests.
- `02-validation/`: simulation drivers, summary CSVs, and validation reports.
- `03-theory/`: technical appendix source.
- `04-application/`: application scripts and aggregate results.
- `05-writing/`: manuscript sources, bibliography, tables, figures, and compiled PDF.
- `review/`: numerical audit, editorial change records, and visual QA.

Raw patient data, replicate caches, and compiled libraries are not included in this revision directory. Application scripts refer to the existing `mapbart-case-study-mm` data directory in the repository. See [distribution notes](review/distribution.md) for packaging details.

## Rebuild

Run from the repository root, with R and a LaTeX installation available:

```sh
export LRC_REVISION_ROOT="$PWD/revision"
R CMD INSTALL revision/01-code/lrcbart
Rscript revision/05-writing/make_tables.R
cd revision/05-writing
latexmk -pdf -recorder -outdir=build -interaction=nonstopmode -halt-on-error main_draft.tex
```

The checked-in summary CSVs support table regeneration without rerunning the simulations. Full simulation and application runs require the dependencies named by their scripts and can be expensive. Historical tuning and development scripts document intermediate analyses; use the full-study drivers and associated reports to interpret the final results.

The manuscript passed numerical and prose-integrity checks and rendered-page inspection before packaging. See [numerical audit](review/numerical_audit.md) and [visual QA](review/visual_qa.md). Remaining scientific questions are documented in the review ledger.
