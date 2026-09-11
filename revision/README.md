# LRC-BART manuscript and code revision

This directory contains the September 10, 2026 manuscript candidate, its implementation, and the summary results used to produce the tables. The original repository files remain available alongside it.

Start with the [paper](05-writing/main_draft.pdf), [cumulative tracked changes](05-writing/main_tracked.pdf), and [explicit f model and ESS clarification](review/f_model_20260911/report.md). The current paper is 80 pages and the cumulative redline is 101 pages. Page 5 displays the BART expansion and prior for f, followed by the discrepancy-leaf mixture in equation (2) and the induced pointwise mixture for f_1. Page 9 explains local ESS as a precision-equivalent reference count. The [46-page author deck](06-slides/main_lrcbart.pdf) includes the rewritten main talk and technical backup; see the [slide rewrite report](review/style_rewrite_20260911/report.md).

The paper now defines a working variance-based ESS, corrects leaf aggregation and the calibration theorem, states implementation overrides, and separates support identification from tree extrapolation. G.6 is a conditional prior-density comparison with a corrected optimizer. CAHB is correctly described as nonparametric, and the application is a descriptive comparison with post-induction ASCT adjustment and an explicit equal-residual-variance assumption.

## Review requested

Yunxuan, please review the corrected scientific definitions and conditional results alongside the producing code, especially the working-ESS interpretation and calibration conventions, the scope of G.6, and the clinical interpretation of the application. The [side-review reconciliation](review/next_revision_20260910/side_review_reconciliation.md) maps all nine findings from the earlier 82-page baseline to the current text.

Two independent agent reviews in each of two rounds found their original major objections resolved, with final minor corrections applied. These are internal AI assessments, not external scientific approval. Author names/order/affiliations are intentionally blank at the user's request for now. The [implementation-impact note](review/next_revision_20260910/implementation_impact.md) records the unused prior-ESS helper limitation and which future changes require new fits or an application reanalysis.

## Contents

- `01-code/`: R/C++ sampler package, comparators, scenarios, and tests.
- `02-validation/`: simulation drivers, summary CSVs, and validation reports.
- `03-theory/`: technical appendix source.
- `04-application/`: application scripts and aggregate results.
- `05-writing/`: manuscript sources, bibliography, tables, figures, and compiled PDF.
- `06-slides/`: slide deck, diagram, bibliography, and source.
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

Validation for this revision: clean and tracked LaTeX rebuilds, complete page overviews and targeted full-size inspections, and bounded deterministic checks against the actual C++ variance calculation. All protected local code/data/result files and numerical table/figure files were preserved. No simulations or application fits were rerun. The original package parsing/install and 102-assertion smoke checks remain earlier evidence, not a new full-model qualification. See [latest presentation QA](review/prior_mixture_order_20260910/visual_qa.md), [its integrity disposition](review/prior_mixture_order_20260910/integrity_disposition.md), [preceding notation QA](review/leaf_prior_clarity_20260910/visual_qa.md), its [integrity disposition](review/leaf_prior_clarity_20260910/integrity_dispositions.md), and the [preceding scientific-round QA](review/next_revision_20260910/visual_qa.md) and [integrity dispositions](review/next_revision_20260910/checks/integrity_dispositions.md).

The [commensurate-center clarification](review/commensurate_center.md) remains preserved; the current induced-mixture slide identifies the all-spike and all-slab components, and the leaf slide explains how each indicator chooses its variance.
