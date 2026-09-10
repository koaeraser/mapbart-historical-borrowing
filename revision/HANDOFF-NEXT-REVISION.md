# Handoff: next LRC-BART revision

Written September 10, 2026. This is the current restart note. It supersedes the old Claude handoff's task status and the status paragraphs in `05-writing/codex_revision_report.md`; retain those files as history. Read the actual current sources before relying on any earlier review verdict.

## Start here

The user is resetting the session and wants another revision round. The last request completed was to make the external-response mean explicit as the center of the commensurate component in both paper and slides. Do not revert to a zero-centered leaf-prior explanation without first showing what the RCT response mean is centered on.

Read in this order:

1. This note and `05-writing/critic_round2.md` / `grader_round2.md`.
2. `05-writing/humanize_logs/codex_20260910/resolution_ledger.md`.
3. Current Section 2.2–2.5 in `05-writing/sections_1_2.tex`, the relevant theory in `03-theory/appendix_new.tex`, and G.6 in `05-writing/web_appendix_G.tex`.
4. `../yunxuan-repo/revision/review/commensurate_center.md` and current slides 5–8.
5. The full current manuscript and slide deck before revising their overall argument.

This is a working candidate for scientific review, not a submission-ready or independently approved paper. Round-2 critic/grader files are internal assessments written by the editing agent, explicitly not independent reviews. The provisional grader result was 33/50, major revision.

## Files and repository

Paths below are relative to the local revision root, even when reading the copy of this note in the Git distribution.

Local revision root:
`/Users/yuanj/Dropbox (Personal)/YuanJi/research/Yunxuan-Zhang/MAP-BART/revision`

| Artifact | Current path |
|---|---|
| Clean paper, 82 pages | `05-writing/build/main_draft.pdf` |
| Paper assembly | `05-writing/main_draft.tex` |
| Main text | `05-writing/sections_1_2.tex`, `sections_3_5.tex` |
| Theory A–F | `03-theory/appendix_new.tex` |
| Appendix G | `05-writing/web_appendix_G.tex` |
| Cumulative colored diff, 90 pages | `05-writing/main_tracked.pdf` and `.tex` |
| Current 25-slide deck | `/Users/yuanj/Bayesics Dropbox/Bayesics Team Folder/Training/JSM-2026-MAP-BART/main_lrcbart.pdf` and `.tex` |
| Diagram source/PDF | Same deck directory, `inserts/lrc_model_diagram.tex` and `.pdf` |
| Git distribution checkout | `../yunxuan-repo` |
| Distribution copies | `../yunxuan-repo/revision/05-writing/` and `06-slides/` |

The revision root is not a Git repository. Commit and push from `../yunxuan-repo`.

- Upstream: `oliviazhang0416/mapbart-historical-borrowing`, default branch `main`.
- Fork: `koaeraser/mapbart-historical-borrowing`; local remote `fork`.
- Branch: `revision/lrc-bart-20260910`.
- PR: https://github.com/oliviazhang0416/mapbart-historical-borrowing/pull/1
- Last verified implementation/upload commit before this handoff: `2018c91261ce1227c3772fbc09915b390c1b3c10`. The handoff itself may add a later documentation commit.
- PR was OPEN, with no reviews or comments, when checked for this note. Recheck after reset.
- Authenticated account `koaeraser` lacks upstream push/reviewer-assignment permission. The user explicitly authorized upload and asking Yunxuan to review. The PR description tags `@oliviazhang0416` and asks her to review. Formal reviewer assignment was rejected by GitHub permissions; no formal reviewer is assigned.
- `note_to_yunxuan.md` was left untouched; do not resend it by assuming the PR request authorizes a separate old message.

## What has been completed

The original handoff's CAHB rows, realized regional ESS panel, control-profile/far-boundary panel, and three approved content moves are integrated. `make_tables.R` generates the tables and `quoted_numbers.txt` from the summary CSVs. Abstract, Sections 1–5 and G received a substantial humanize-prose pass. G.6 was retained pending scientific review.

Keep these distinctions:

- CAHB's reported precision share is not LRC's posterior compatibility probability or realized ESS. Its main Table 2 diagnostics were suppressed where definitions differ; the appendix reports its own diagnostic definitions.
- CAHB uses R=200 paired replicates; the primary full-study rows use R=500. The transferred implementation is ours. Avoid claiming a general method ranking or a demonstrated dimensionality limit from this comparison.
- Table 2 surface RMSE uses all RCT profiles. Control-profile and far-boundary summaries use their explicitly stated profile sets.
- Regional ESS values use R=200 and regional standardization; they are not additive across regions. Sc1: 32.7/31.9; Sc4 d1: 15.1/12.7; Sc4 d2: 9.5/0.44; Sc5 d2: about 0.2/0.2. Verify against `regional_ess_summary.csv` when quoting.
- Removing boundary profiles does not eliminate the Sc4 surface loss: d2, far compatible region, LRC 0.747 versus BART-PP 0.723, paired difference 0.024 with SE 0.005.
- Application calibration concerns the control log-time mean, not RMST itself.
- The sampler truncates the spike-scale prior. Student-t marginal and zero-limit/convexity claims must retain the stated untruncated-law qualifications.

### Latest correction: external-mean centering

Paper Section 2.2 now defines f as the external-control response mean and f_1=f+g as the RCT-control response mean. Before unchanged equation (2), it displays the conditional all-spike distribution

`f_1(x) | f, T_g, {z_h,ell_h(x)=1 for all h}, tau_0^2 ~ N(f(x), H_g tau_0^2)`.

It then introduces g=f_1-f, explaining the zero-centered leaf prior. The general conditional variance is the sum over reached leaves of `z*tau_0^2+(1-z)*tau_1^2`. Both components retain the external mean as the center; the slab permits larger departures. The statements concern the latent mean, with observation residual variance separate. Conditions on partitions, states and scales are essential. Do not equate one leaf's variance tau_0^2 with the variance of the entire H_g-tree discrepancy.

Slides: 5 defines f and g; 6 gives the external-centered prior and its discrepancy parameterization; 7 is the tree diagram using the old slide's visual style; 8 explains the leaf mixture. Leaf indicators vary; scales and w are shared. The old `main.tex` / `main.pdf` MAP-BART deck is untouched. No sampler or result changes accompanied this clarification.

## Priorities for the next revision

1. **ESS aggregation across leaves.** Section 2.4 uses a Binomial(H_g,w) expression for a standardized estimand, although the implementation assigns independent indicators to all leaves. Establish the correct fixed-partition variance and averaging over leaf states when the estimand spans multiple leaves. A pointwise profile reaches one leaf per tree; a standardized mean can span many. Do not confuse the two. Trace the calculation to `01-code/lrcbart/R/lrcbart.R`, the calibration functions and Appendix D.
2. **ELIR interpretation.** The paper still states that a scalar prior with variance V has expected ELIR ESS sigma_1^2/V. Determine when this identity is exact and whether it applies to the nonnormal BART/mixed marginal prior used here. Distinguish conditional Gaussian, inverse-variance approximation and exact marginal information. If a change affects calibration, identify the required downstream code/results work before making cosmetic substitutions.
3. **Support and identification.** Section 2.2 still says g is a prior draw where the RCT has no support. Check sharing through trees and learned hyperparameters, and distinguish fixed-partition conditional claims from full-posterior identification. Existing qualifications elsewhere do not resolve this sentence automatically.
4. **G.6 prior-cost argument and its main-text claims.** The optimizer/threshold argument and inference about indicator behavior remain unresolved. Section 2.2 still claims no clinically sized shift flips an indicator; the slide “Why a separate discrepancy ensemble” makes similarly strong statements. Audit the mathematics and narrow claims to what is established. The protected earlier judgment is not scientific approval.
5. **Full slide/paper reconciliation.** The latest slide change was local, not a full synchronization. Recheck old “one residual SD,” “sample-size limit,” surface/map/ESS wording and comparator caveats. Do not assume the new deck has inherited every manuscript correction.
6. **Length, author field and prose.** `\author{}` remains empty. Current texcount text totals: 4,040 (Sections 1–2, including abstract) + 3,774 (Sections 3–5) = 7,814, approximately 7,586 excluding the unchanged 228-word abstract. This is above the approximate 7,000-word target. Recount after changes. Fix the residual comma before “The sampler” in Section 2.2 when editing that paragraph. Preserve needed caveats while cutting.

Resolve scientific issues before another broad prose pass. Verify the actual source rather than treating prior completion reports as evidence that all problematic claims were removed. If independent review is run, label it accurately and keep it distinct from the editing agent's assessment.

## User preferences and constraints

Use `/Users/yuanj/.codex/skills/humanize-prose/SKILL.md`, including its house-voice and slide-mode references. The user dislikes Claude's prose and wants natural, precise authorial exposition. No production/agent/reviewer narration in the manuscript. Explain what is being borrowed and at which level. Avoid label hooks, metaphorical claims and unexplained zero-centered priors.

Snapshot before editing; preserve the original and the current reviewed candidate. Preserve numerical tables/results, code and Appendix A–F unless a specific scientific revision requires a documented change. The earlier authorized Appendix D implementation remark is already included. Do not silently rerun, replace or omit results. The open scientific issues are not a reason to rewrite equations into unverified certainty.

## Build, validate and upload

From local `05-writing`:

```sh
latexmk -pdf -recorder -outdir=build -interaction=nonstopmode -halt-on-error main_draft.tex
latexmk -pdf -recorder -outdir=build -interaction=nonstopmode -halt-on-error main_tracked.tex
```

The tracked source is generated, not an input to the clean paper. After a tracked rebuild, copy `build/main_tracked.pdf` to `main_tracked.pdf` for the stable delivered path. From the deck directory, use the same latexmk command without `-outdir=build`, targeting `main_lrcbart.tex`. A stale auxiliary state once made latexmk invoke BibTeX before citations existed; a successful direct pdflatex pass followed by latexmk fixed it. The deck retains existing Metropolis/pdfLaTeX and font-substitution warnings; there are no overfull boxes or unresolved references.

Last validation: clean paper and tracked paper have no warnings or overfull boxes. Changed pages/slides were checked at full size; all 82 clean-paper pages and 25 slides were inspected in overview sheets. The tracked model pages were checked separately. Initial prose-only comparisons against `post_technical/` passed strict integrity for all three manuscript files. The latest centering pass adds explanatory math and one reference to equation (1); its checker exception and audits are recorded under `../yunxuan-repo/revision/review/commensurate_*`.

At initial repository packaging, all 53 R files parsed, the package compiled/installed, and the single-leaf/two-leaf tests passed 102 assertions. Full simulations were not rerun. Do not describe those tests as full scientific validation.

The distribution code uses `LRC_REVISION_ROOT` (default `revision` from repository root) in eight scripts in place of local paths; application data paths refer to the repository's existing data directory. See `revision/review/distribution.md`. Do not overwrite these portable path edits when syncing local files. Large replicate caches, local patient-data copies and compiled libraries are intentionally excluded.

Sync changed manuscript sources/PDFs to `../yunxuan-repo/revision/05-writing/`; slides, bibliography and required inserts to `revision/06-slides/`. Inspect the diff, commit, push to `fork revision/lrc-bart-20260910`, and verify the PR's current head. No need to ask again for the previously authorized upload or review request.

## Baselines and tracked-diff lineage

- Original pre-Codex paper baseline: `05-writing/humanize_logs/codex_20260910/` (main source/PDF, sections, original tables and theory).
- Prose-validation baseline: `.../codex_20260910/post_technical/`; this was a technical-edit checkpoint, not scientific certification.
- Before the centering correction: `05-writing/humanize_logs/commensurate_center_20260910_092645/`.
- Restart baseline captured for this handoff: `05-writing/humanize_logs/handoff_20260910_093747`. It includes current manuscript/slides, tables, theory, prior handoff and `SHA256SUMS.txt`. Do not edit this snapshot.

The cumulative colored diff compares the original pre-Codex manuscript with the current paper. It has blue underlined additions and red struck deletions. Tables 1–2 appear as full red “before” and blue “after” versions; two new appendix tables are blue. Other tables show current formatting. The bibliography is current and is not individually tracked.

Full `latexdiff --flatten` breaks changed table/longtable structures. The successful process expands the main text/theory/appendix inputs while replacing table inputs with stable placeholders, runs `latexdiff --math-markup=whole --disable-citation-markup`, then restores table content with explicit before/after treatment. It also repairs diff markup around the formatting-only emergency-stretch assignment, hides hyperlink borders, and uses a local sloppypar for the first changed model paragraph. The handoff snapshot contains copies of `tracked_prepare.py` and `tracked_finish.py` as generation aids; they currently reference `/tmp/lrc-tracked` and require path adaptation and those final formatting fixes. Temporary directories are not durable baselines. The complete delivered `main_tracked.tex` builds directly using the current bibliography and figures.
