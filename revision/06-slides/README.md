# LRC-BART presentations

`main_lrcbart.pdf` is the rewritten 46-page author presentation. The 33-page main talk is followed by one references page and 12 pages of technical backup. The visual and writing style follows Yuan Ji's 2022 ENAR BaySize deck. The leaf mixture is on slide 9, the explicit indicator formulation on slide 10, the induced pointwise mixture for $f_1$ on slide 11, and the discrepancy-tree diagram on slide 8.

Build the author deck here with `latexmk -pdf -recorder -interaction=nonstopmode -halt-on-error main_lrcbart.tex`. Existing pdfLaTeX font substitutions are retained.

`discussants/` contains two editable PowerPoint decks: an industry discussion and a hypothetical FDA-perspective discussion. Each has six core slides for 5–7 minutes, two backup slides, speaker notes and sources. The FDA deck represents a hypothetical discussant view, not FDA policy or endorsement. It distinguishes the January 2026 Bayesian draft and February 2023 external-control draft from final ICH E9(R1) guidance.

The substantial style-rewrite audit is in `../review/style_rewrite_20260911/`. The earlier prose-only pass is retained in `../review/slides_voice_20260911/`.

The deck defines BART-PP as pooled-control BART with a source indicator. It presents RMST results across all five main scenarios with a concise median-ratio caveat. The application section reports the primary analysis without discussion of prior coding or debugging history. See `../review/survival_balance_20260911/`.
