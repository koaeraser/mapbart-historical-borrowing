# LRC-BART next revision: completed scientific and editorial round

September 10, 2026. The scientific-round report below records the candidate before the subsequent notation follow-up. Read [the prior-mixture presentation](prior_mixture_order.md) for current artifacts: 79-page clean paper, 101-page cumulative redline, 25 slides, and 6,932 main-text words excluding the abstract. The archived scientific reviews and their scope remain unchanged.

The manuscript now defines the implemented variance-based ESS, uses independent leaf-state aggregation, corrects the calibration theorem and its code exceptions, separates support from extrapolation, and replaces G.6 with a scoped constrained prior-density calculation. CAHB is correctly positioned as nonparametric. The single-arm application is descriptive with post-induction ASCT adjustment and an explicit equal-residual-variance prediction assumption; nonlinear survival means are no longer claimed invariant. The external-control mean remains the center of the commensurate component. Slide 6 now labels the all-spike conditioning and explains where its indicators enter.

The main text contains 7,220 texcount text words including the 227-word abstract: 6,993 excluding the abstract, versus approximately 7,586 at restart. This count excludes separate table files, mathematical-expression counts, headings and captions. The clean paper is 79 pages; the cumulative tracked paper is 99 pages; the deck retains 25 slides.

Two method/applied agent reviews ran independently in each of two rounds. Each second-round reviewer marked all five of its original major comments resolved, reported zero major comments and requested minor corrections, which were applied in the ledger. The user-supplied side review covered the pre-run 82-page baseline; all nine major findings were reconciled against the candidate, with residual wording/assumption fixes applied. These are internal AI assessments, not external peer review or submission approval.

Validation:

- Original and restart baselines are immutable snapshots. All 304,665 protected code, validation and application files are unchanged; numerical table and figure files match the restart baseline.
- Bounded deterministic checks compile and exercise the actual C++ leaf-variance function and verify the corrected algebra. No simulation or application results were rerun.
- Clean, cumulative and per-round tracked builds succeed without LaTeX warnings, unresolved references/citations or overfull boxes. The deck has no unresolved references/citations or overfull boxes, with existing Metropolis/font substitution warnings retained.
- Complete rendered-page overviews and targeted full-page/slide checks were performed. Redline old-table color inheritance was corrected after visual inspection.
- Prose-only checkpoint checks pass; verified scientific deltas are explicitly mapped in checks/integrity_dispositions.md rather than represented as a blanket prose-only pass.

Read `side_review_reconciliation.md`, `implementation_impact.md`, the two round ledgers, `source_verification.md` and `visual_qa.md` for scope and evidence. The unused exported lrc_ess_prior approximation and outdated roxygen text remain documented code debt. Exact marginal ELIR, changed calibration rules or a changed application adjustment set would require new empirical analyses. The former MAP-AFT-BART sampler remains unavailable; earlier tail differences have no isolated causal attribution.

The user has requested that author names, order and affiliations remain blank for now; this is intentional and no clarification is pending. Independent clinical qualification, full sampler convergence and confirmatory operating characteristics remain outside this work.

## Artifact lineage

- Current local paper: `05-writing/build/main_draft.pdf` assembled by `05-writing/main_draft.tex`.
- Current local cumulative redline: `05-writing/main_tracked.tex` and `.pdf`, against the original pre-Codex snapshot `05-writing/humanize_logs/codex_20260910/`.
- Current deck: `main_lrcbart.tex` and `.pdf` in the JSM-2026-MAP-BART training folder.
- Run: `05-writing/critical_revisor_logs/run_20260910_094551/`.
- Restart baseline and SHA manifest: `round_0/`.
- Round-1/round-2 clean sources, PDFs, reviews, ledgers and redlines: corresponding round folders. Build archived clean sources from their `clean/05-writing/` directory.
- `build_tracked.py` records the table-aware custom redline workflow used because full latexdiff flattening breaks the supplied table/longtable structures. Bibliography and typography are current and untracked. Main Tables 1–2 have full before/after treatment; new cumulative appendix tables are blue.
- Git distribution: `../yunxuan-repo/revision/`, existing branch `revision/lrc-bart-20260910`, PR https://github.com/oliviazhang0416/mapbart-historical-borrowing/pull/1 . The upload authorization in the prior handoff persists; portable R paths and note_to_yunxuan.md are preserved. No separate message is sent.
