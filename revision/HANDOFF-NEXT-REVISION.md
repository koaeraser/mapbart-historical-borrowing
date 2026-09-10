# Handoff after the next LRC-BART revision

September 10, 2026. This is the current restart note. The previous handoff is preserved in `05-writing/critical_revisor_logs/run_20260910_094551/round_0/`. Earlier editing reports and their verdicts remain historical assessments of earlier versions.

## Current state

The next scientific revision is complete, with two independent method/applied agent reviews in each of two rounds, followed by integration of the user's `FULL-PAPER-REVIEW.md`. Each round-2 reviewer marked all five of its prior major comments resolved and requested minor corrections, which were applied. The side review assessed the pre-run 82-page baseline; its nine major findings are reconciled in the current report. These are internal AI assessments, not external peer review or scientific/clinical approval.

Read first:

1. `05-writing/leaf_prior_clarity.md`, then `05-writing/next_revision_report.md`.
2. `05-writing/side_review_reconciliation.md` and the unchanged `FULL-PAPER-REVIEW.md`.
3. `05-writing/critical_revisor_logs/run_20260910_094551/implementation_impact.md`.
4. Both round ledgers and final sources in that run.

Current artifacts:

| Artifact | Path |
|---|---|
| Clean paper, 79 pages | `05-writing/build/main_draft.pdf` |
| Main source | `05-writing/main_draft.tex`, `sections_1_2.tex`, `sections_3_5.tex` |
| Theory A–F | `03-theory/appendix_new.tex` |
| Appendix G | `05-writing/web_appendix_G.tex` |
| Cumulative tracked changes, 101 pages | `05-writing/main_tracked.tex`, `main_tracked.pdf` |
| 25-slide deck | `/Users/yuanj/Bayesics Dropbox/Bayesics Team Folder/Training/JSM-2026-MAP-BART/main_lrcbart.tex` and `.pdf` |
| Run, baselines, reviews, checks, renders | `05-writing/critical_revisor_logs/run_20260910_094551/` |
| Git distribution | `../yunxuan-repo/revision/` |

Main text: 7,188 texcount text words including the 227-word abstract; 6,961 excluding it. Separate tables, mathematical-expression counts, headings and captions are excluded.

The user's subsequent notation follow-up is complete: equation (2) defines the leaf prior with variance tau_{1-z_hl}^2, theta_hl and the reached-leaf index are explicit, and all-slab/mixed f_1 priors are displayed. Appendix A/B and slide 8 use the same indexing. See `05-writing/leaf_prior_clarity.md` and its separate follow-up run; earlier review snapshots remain unchanged.

## Scientific decisions to preserve

- The commensurate component centers the RCT control response mean f_1 at the external-control mean f. The zero-centered discrepancy is g=f_1-f. The displayed variance H_g tau_0^2 is conditional on every reached leaf indicator being 1. The user asked where z enters; slide 6 now explicitly labels the all-spike case and explains that each z selects tau_0^2 versus tau_1^2. Do not describe that display as the general prior marginalized over z.
- ESS is a defined variance-based working precision summary. Only the normal reference prior yields the exact marginal ELIR variance ratio. No predictive-consistency claim is made for the BART mixture.
- A standardized mean can span multiple independently labelled leaves; V_g is the sum of squared target leaf proportions times state-specific variances. Binomial(H_g,w) is restricted to one positive-weight leaf per tree. Existing realized/regional results already use the exact weighted-leaf function.
- Calibration retains the implemented mean partition coefficient before inversion, untruncated scale-law integration, and a truncated fitting prior. The rescaled block-law theorem concerns exact integration and grid selection under explicit assumptions. The ceiling cap and upper-grid fallback are separate code branches. A change in these conventions would require new fits.
- The ceiling includes finite BART regularization. Its external-sample-size equality and covariate-shift reduction are a flat-prior, fixed-one-tree reference, not universal finite-prior claims.
- Population support identification, shared-tree extrapolation and conditional empty-leaf priors are distinct. Only the globally single-arm g process is likelihood-independent everywhere.
- G.6 retains all spike contributions in the constrained optimizer, uses integer slab counts and is a conditional prior-density comparison. It implies neither clinical indicator-switching thresholds nor inability of shared-leaf models to borrow locally.
- CAHB has nonparametric kernel outcome/precision estimates. Its transferred R=200 results and diagnostic definitions retain the documented caveats versus the R=500 primary study.
- Gaussian linear means can remain centered while nonlinear survival-ratio means move under discrepancy uncertainty. The application w=0.9 PFS mean changes from 0.97 to 1.00. Empirical near-zero simulation bias is not called theoretical unbiasedness.
- The application adjusts for five baseline covariates plus ASCT after induction. It is a descriptive model-based association conditional on recorded profiles, not a total causal treatment effect. Equal control residual variances are an explicit assumption for counterfactual survival prediction as well as calibration. Changing the adjustment set or variance assumption requires application reanalysis.

## Preserved results and open boundaries

All 304,665 protected files in 01-code, 02-validation and 04-application were verified unchanged. Numerical table/figure files and the diagram were preserved. No simulation or application results were rerun. Bounded deterministic checks exercise the actual C++ weighted variance and revised mathematics; they do not certify the full sampler or operating characteristics.

The unused exported `lrc_ess_prior` helper still uses the tree-binomial approximation; old roxygen summaries are documentation debt. It is not used to produce the reported results. See implementation_impact.md before changing code or interpreting this helper as corrected.

The original MAP-AFT-BART sampler remains unavailable, so earlier tail differences cannot be assigned to a particular prior, coding or sampler change. Clinical qualification, full-model convergence and confirmatory operating-characteristic work remain outside this revision. The user explicitly requested that author names, order and affiliations remain blank for now. This is intentional and is not a pending clarification. Do not infer author metadata from slide credits or ask again unless the user reopens authorship.

Dense simulation tables, their continuation/regional panel and float whitespace remain optional presentation work. Hyperlink borders were hidden. No rows/figures were dropped or redesigned during this round.

## Build and artifact lineage

From 05-writing:

```sh
latexmk -pdf -recorder -outdir=build -interaction=nonstopmode -halt-on-error main_draft.tex
latexmk -pdf -recorder -outdir=build -interaction=nonstopmode -halt-on-error main_tracked.tex
cp build/main_tracked.pdf main_tracked.pdf
```

Compile `main_lrcbart.tex` in the deck folder with latexmk -pdf. Current clean/redline builds have no LaTeX warnings, unresolved references/citations or overfull boxes. Slides retain Metropolis/font substitution warnings, with no overfull boxes or unresolved references. The preceding scientific run contains complete overview inspection; the separate leaf-prior follow-up contains targeted final-page inspection in its visual_qa.md.

The cumulative redline compares the original pre-Codex baseline in `05-writing/humanize_logs/codex_20260910/` with the current paper. `run_20260910_094551/build_tracked.py` expands text/theory inputs using table placeholders, runs latexdiff with whole-math markup, and restores full red-before/blue-after main Tables 1–2 plus blue new cumulative appendix tables. Old captions retain black text with Before revision labels. Bibliography and typography are current and untracked. Round-specific redlines and clean checkpoints are archived separately. Never overwrite immutable baselines.

## Repository and retained authorization

The revision root is not Git. Work from `../yunxuan-repo`, branch `revision/lrc-bart-20260910`, fork remote `koaeraser/mapbart-historical-borrowing`, upstream `oliviazhang0416/mapbart-historical-borrowing`. Existing PR: https://github.com/oliviazhang0416/mapbart-historical-borrowing/pull/1 . Recheck live state before asserting head/review status.

The user previously authorized upload and asking Yunxuan to review. Continue syncing the corrected manuscript/slides and curated review record, committing and pushing the same fork branch without asking again. Preserve the portable LRC_REVISION_ROOT edits in distribution R scripts. Do not copy local caches, compiled libraries, patient data or the large protected-file manifest into Git. Do not send a separate note to Yunxuan; note_to_yunxuan.md remains untouched. Formal reviewer assignment previously failed because the authenticated account lacks upstream permission.
