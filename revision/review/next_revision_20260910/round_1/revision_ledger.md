# Round 1 revision ledger

These are dispositions by the integrating editor following two independent AI reviews. They are not external scientific acceptance.

| Finding | Disposition | Location and change | Evidence |
|---|---|---|---|
| R1 M1; R2 m2 | Applied | Sec2.4, D, abstract/intro/deck: working precision, normal reference, nonnormal marginal distinction, uncensored log-time units | ELIR primary paper eq7; code; deterministic mixture check |
| R1 M2 | Applied in manuscript | Sec2.4,D: squared leaf proportions and general independent-state average; realized formula matches production code | R and C++ producing functions; compiled direct check |
| R1 M2 helper | Source-flagged follow-up | Unused exported lrc_ess_prior still computes tree-state shortcut; no numerical outputs depend on it. Documented in code-impact note; numerical code preserved per handoff | Call-site search |
| R1 M3 | Applied | D and main theorem: rescaled block law, finite-m moments/threshold conditions, no unsupported CLT, explicit cap/fallback, flat-prior support conditions | Self-contained proof, code cap branch, anchor counterexample |
| R1 M4; R2 M3 | Applied | Sec2.2,E: population identification separated from shared-leaf/hyperparameter extrapolation; conditional single-arm variance | Model likelihood and prior factorization |
| R1 M5; R2 M1 | Applied | G.6, intro/model/deck: constrained integer prior-density cost, retain spike contributions, remove clinical threshold/posterior guarantee | Lagrange multipliers; deterministic checks |
| R2 M2 | Applied | Related work/G.0: CAHB is kernel-based nonparametric historical control borrowing | Primary CAHB pp5341/5343 |
| R2 M4 | Applied | Sec4,G.5/deck: five baseline covariates plus ASCT after induction; descriptive model-based comparison | prep_data.R and cohort script |
| R1 m3; R2 M5 | Applied | Single-arm/application/G.7/deck: linear mean centering does not preserve nonlinear ratio means | Preserved CSV rows in checks/quoted_corrections.json |
| R1 m1 | Applied | Sec2.4,D, app caption/deck: pointwise plug-in mean scale distinguished from prior-averaged reciprocal | lrc_ess_map and pointwise_ess.R |
| R1 m2 | Applied | Main leaf theorems: interior w and slab comparator | Appendix B conditional formulas |
| R1 m4 | Applied | Appendix A kernel-invariance wording | MH invariant-kernel construction |
| R2 m1 | Applied | Full deck scientific reconciliation; original centering diagram retained | Deck audit; CSV panels |
| R2 m3,m4,m6 | Applied | Coding extrapolation, cohort masks, survival time notation, C-BART description, log-RMST diagnostics | Producing R scripts |
| R1 m5; R2 m5 | Partially applied | Main text shortened; D self-contained; obsolete D cross-references removed. Historical A/C/F prose outside scientific dependencies remains protected | Handoff A-F preservation constraint; prose checker |
| Source flag: old tails | Applied by narrowing | G.5/deck: observed comparison without causal attribution to old prior | Old sampler unavailable |
| Source flag: author field | Pending user metadata | Blank field preserved; names/order/affiliations requested asynchronously | Original source also blank |

Post-technical to shortened-main hard integrity checks pass numbers, citations, labels, sections, numbered equations, protected blocks, graphics, references and structure. Inline-math audit removes repeated mentions of M_f, M_g, g and T_i, whose definitions remain in equations/prose. Number-word changes remove repeated narrative counts; all design counts remain stated. Tables, figures and data are unchanged.
