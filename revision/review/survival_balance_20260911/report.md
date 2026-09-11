# Comparator clarification and survival-slide balance

BART-PP is pooled-control BART with source membership as a predictor, not a power prior. Its definition follows manuscript Section 3 and Web Appendix G.2.

Only slides 14 and 18 changed. Slide 18 now shows LRC-BART (target 100) and BART-PP RMST bias, RMSE and coverage for all five main survival scenarios. Values were extracted directly from the current manuscript table. The slide leads with gains under agreement and retains a concise caveat about median-ratio performance under unmeasured confounding and 90% regional RMST coverage.

The prior 25-slide deck is preserved in pre_humanize.tex/pdf. All other source frames are byte-identical, and PDF text differs only on slides 14 and 18. The slide-mode integrity checker passes with the explicitly authorized table and explanatory changes. Numeric, inline-math and table deltas are confined to the authorized survival display and were checked against the producing manuscript table. Citations, figures, labels, sections and environment structure are unchanged.

Clean build: 25 slides, no overfull boxes or unresolved citations/references. Existing font/theme substitutions remain. Both changed slides were rendered for visual inspection.
