# LRC-BART author-deck rewrite

The author deck was rewritten from the immutable 25-slide baseline in `pre_style_rewrite.tex`, using Yuan Ji's 2022 ENAR BaySize deck as the primary style reference. The final PDF contains a 33-page main talk, one references page and 12 pages of technical backup. Short section breaks, ordinary bullets, larger equations and figures, and brief interpretation replace the earlier dense block layout.

The rewrite preserves the manuscript's distinctions among compatibility probability, working ESS, pointwise ESS and realized regional ESS. It also preserves the conditional scope of the theoretical results, the descriptive interpretation of the EloKRd comparison, and the inability of single-arm data to identify a control discrepancy. Comparator labels and the Gaussian and survival scenario grids were checked against the paper sources.

At the presenter's request, the application no longer discusses the discovered coding problem, old results, the former long tail or the unavailable sampler. Slide 29 now presents the application analysis. Slide 31 explains that red dots are PFS, blue dots are OS, each dot evaluates the external-information ESS at one EloKRd profile, and categorical jitter has no statistical meaning.

The PDF was rebuilt with `latexmk -pdf -recorder -interaction=nonstopmode -halt-on-error main_lrcbart.tex`. It has 46 pages at 16:9 size, with no overfull boxes, unresolved references, unresolved citations or LaTeX errors. All pages were rendered and inspected; slides 29 through 31 were re-inspected after the presentation follow-up. `source.diff` is the cumulative source review artifact.
