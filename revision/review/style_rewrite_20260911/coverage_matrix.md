# LRC-BART deck style rewrite coverage

Reference: `presentation/invited/ENAR-2022/BaySize Talk/Bayesize-v1.tex` and its 34-page PDF.

Authorized scope: substantial rewrite of the LaTeX slides to reduce density and match Yuan Ji's pre-2024 presentation style. Frame insertion, deletion, splitting, reordering, and movement of technical material to backup are allowed. Scientific claims, estimates, equations, tables, figures, citations, and caveats remain fixed except for source-verified clarifications listed below.

| Old slide | Narrative job | Decision |
|---:|---|---|
| 1 | Title | revise to UChicago 2022 visual style |
| 2 | Overview | replace with short talk map |
| 3 | Motivation for external controls | split gain and local mismatch |
| 4 | Existing borrowing priors | simplify to a three-row comparison |
| 5 | Outcome model | retain equation; shorten definitions |
| 6 | Leaf mixture | split mixture from interpretation of z |
| 7 | Induced prior for trial control mean | retain equation; shorten interpretation |
| 8 | Model diagram | retain as a visual explanation |
| 9 | Separate ensembles | shorten; move density qualification to backup |
| 10 | Tree updates | move detailed factors to backup |
| 11 | Prior ESS | split working ESS, ceiling/rule, and realized diagnostic |
| 12 | Theory | move detailed theorem conditions to backup; retain concise main-slide summary |
| 13 | Borrowing map | retain definition and interpretation; shorten conditions |
| 14 | Simulation design | split settings from scenario definitions |
| 15 | Gaussian results | simplify main comparison; retain full table in backup |
| 16 | Regional borrowing map | enlarge figure; reduce prose |
| 17 | Regional conflict | split performance loss and regional ESS |
| 18 | Survival results | simplify main table; retain full metrics in backup |
| 19 | Single-arm results | simplify to identification limit and key numbers; retain details in backup |
| 20 | EloKRd and UCMM data | retain paired source summary; present the application specification without coding history |
| 21 | EloKRd results | retain the primary table and interpretation; omit historical sampler and debugging details |
| 22 | Local ESS map | add a visible endpoint key and explain profiles, variance-equivalent information and jitter |
| 23 | Main findings | reduce to three takeaways |
| 24 | Limitations | reduce to four audience-facing limits |
| 25 | References | retain bibliography |

## Source-verified clarifications

1. State the Pr-rule direction explicitly: working ESS is at or below the target in at least 95% of blocks; a target at or above the whole-chain ceiling triggers the smallest-scale override.
2. Separate Gaussian and survival regional/global shift grids.
3. Describe `H_f=50, H_g=5` as simulation settings; the application uses `H_f=10, H_g=5`.

## Semantic locks

- Compatibility probability, prior ESS, pointwise ESS, and realized regional ESS remain distinct.
- The borrowing-map guarantees remain conditional on fixed partitions and fixed `f`.
- Exact collapsed tree-move factors do not establish finite-chain mixing.
- The EloKRd comparison remains descriptive because ASCT is measured after induction.
- Single-arm data cannot identify a control discrepancy without internal controls.

## Presentation-only follow-up

- Slide 29 no longer reports the discovered ethnicity-coding problem or before-and-after ESS ceilings. It now introduces the application estimand, tree counts, ESS targets and posterior summaries.
- The technical-backup slide on the earlier EloKRd analysis, former long tail and unavailable sampler was removed.
- Slide 31 now identifies red as PFS and blue as OS. It states that each dot is the local ESS for one trial profile, that higher values indicate more external-control information, and that categorical jitter only separates overlapping points.
