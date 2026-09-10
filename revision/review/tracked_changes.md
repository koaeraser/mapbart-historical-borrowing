# Colored tracked manuscript

The 90-page `05-writing/main_tracked.pdf` compares the preserved pre-revision manuscript with the September 10 candidate. Blue underlined prose denotes additions; red struck prose denotes deletions. Changed equations are tracked as whole equations.

Tables 1 and 2 appear as full before (red) and after (blue) versions because their table structure changed. The two new appendix tables are blue. Other tables use the current layout and unchanged numerical results. Formatting changes and the bibliography are not individually tracked. Relocations may appear as deletions and additions.

The self-contained LaTeX source embeds the compared manuscript text and tables and uses the existing `figures/` and `merged.bib`. From `revision/05-writing`, run:

```sh
latexmk -pdf -outdir=build -interaction=nonstopmode -halt-on-error main_tracked.tex
```

Generated with latexdiff, with explicit table handling and a repaired formatting-only emergency-stretch command. The final build has no warnings, undefined references, or overfull boxes. All 90 pages were inspected in rendered overview sheets; the before/after tables were also checked at larger scale.

Baseline SHA-256 hashes (preserved local `humanize_logs/codex_20260910` snapshot):

- `main_draft.tex`: `f578e10e2c6b7e020b4923509bd6b3e8edb9354bf0c5bc32bc8d14231d1dd372`
- `sections_1_2.tex`: `d99ee98bef83521a0cd5173b21a9667e2854b61bf7508cad6730b070b88e8703`
- `sections_3_5.tex`: `c48b6f3a529cfedcceca602dd8353db25cef95f237ebb7848350f1c5dfcd1f2f`
- `web_appendix_G.tex`: `d57a0fa3a51ce3957680127e50198c8b120af6454ca2a87edc831b69768e6d2a`
- `appendix_new.tex`: `933a23c63e2cd9ef0a080bfe7c46e97d6a893024c43ca8ba4d470216534b3525`
