# Summary

Independent AI methodologist review of the September 10 restart manuscript assembled by `05-writing/main_draft.tex`, both main-text section files, all input tables, `03-theory/appendix_new.tex`, and `05-writing/web_appendix_G.tex`. I read the complete text before formulating the review and did not read the contextual review. This is a scientific manuscript review, not external human approval or a complete implementation/convergence audit. No manuscript, result, or model code was edited. Algebra checks are in `reviewer1_algebra_checks.txt` beside this review.

The model and conditional leaf results are substantially clearer than the calibration and design-justification claims. Five actionable issues block the current mathematical exposition. The main numerical results can largely be preserved by defining the implemented precision summary honestly, aligning the stated realized summary with its producing code, and narrowing theorems and interpretation. Replacing the actual calibration algorithm, in contrast, would require downstream refits.

# What the paper does well

- The external response mean is explicit as the center of the commensurate component; the zero-centered discrepancy leaves then have an intelligible role.
- Appendix A gives a valid marginal-then-conditional kernel argument when its stipulated MH ratios and update order are used. The code inspected redraws a tree's leaves immediately after its collapsed move (`lrcbart.cpp:343-416`) and then updates hyperparameters after the ensembles (`:674-681`). This supports the specified algorithmic ordering, without proving that all proposal details are correct or that finite chains have mixed.
- The leaf influence bound clearly separates the additional spike contribution from ordinary slab shrinkage. The fixed-partition ensemble consistency proof states representability and avoids identifying mixture indicators with the substantive compatibility map.
- The simulation discussion admits moderate-conflict bias, imperfect map separation, surface loss away from the boundary, and the transfer-specific limits of CAHB. These are valuable and should survive shortening.

# Major comments

## M1. Define the information summary being calibrated; it is not general marginal ELIR

**Location.** `sections_1_2.tex:122-133`; abstract and introduction; Appendix D, especially `appendix_new.tex:382-387,436-443`.

**Issue and consequence.** The assertion that a scalar prior of variance V has expected ELIR equal to sigma1^2/V is false for a general prior. The primary paper distinguishes variance/precision ratios in Section 2.2 from ELIR in Section 2.4, equation (7): ELIR averages the negative second derivative of the log prior divided by unit Fisher information. Equality with inverse variance holds for a normal prior and normal mean likelihood with fixed residual variance. Its Student-t example explicitly differs. See [Neuenschwander et al., primary preprint](https://arxiv.org/pdf/1907.04185), Sections 2.2, 2.4 and 2.5.1. For a t3 prior of unit scale and a unit-variance normal observation, inverse marginal variance is 1/3, exact ELIR is 2/3, and averaging Gaussian conditional precisions gives 1.

The implemented target averages reciprocal Gaussian working variances over the spike scale, with a BART posterior represented by its variance and random partitions represented by c_g. Thus it is neither inverse variance of the complete marginal prior nor its exact marginal ELIR. The same distinction applies to the ceiling. Moreover, calling `sigma1^2` a unit information quantity for a standardized regression mean implicitly chooses a normal location reference experiment; it is not the efficient observed-data information for that functional in a general regression or censored survival model.

**Concrete fix.** Retain the current target as a clearly defined working precision ESS, motivated by conditional Gaussian matching. State the normal reference experiment and distinguish the average conditional precision from exact marginal ELIR. Keep the verified ELIR citation as background, without claiming predictive consistency or exact ELIR for this BART mixture. Label the ceiling as the ceiling of this summary, including the regularization in the RWD fit. For survival it is in uncensored log-time normal-reference units and does not calibrate RMST information. Rewrite the abstract/introduction claims consistently.

**Downstream work.** Purely identifying the currently computed summary requires no reruns. Choosing exact marginal ELIR, or the inverse of the complete marginal variance, changes the method: develop and verify the density/information calculation, recalibrate all scales, refit affected LRC-BART/C-BART variants and application ladders, and regenerate their outputs before changing numerical claims.

## M2. Standardized estimands require independent leaf-state aggregation, and realized ESS already uses that aggregation

**Location.** `sections_1_2.tex:122-141`; `appendix_new.tex:387,402,425,443`. Producing code: `01-code/lrcbart/R/lrcbart.R:466-508`, `01-code/lrcbart/src/lrcbart.cpp:754-786`; `01-code/R/fit_lrcbart.R:110-116`; `02-validation/sc4_boundary/pointwise_ess.R:95-113`.

**Issue and consequence.** The prior has one independent indicator per leaf, while the displayed Binomial(H_g,w) sum effectively gives each entire tree one state. Let a_hl be the fraction of standardization profiles in leaf l of tree h. Conditional on partitions, states and scales, the correct prior discrepancy variance is

`V_g(T,z,tau) = sum_h sum_l a_hl^2 [z_hl tau0^2 + (1-z_hl) tau1^2]`.

The all-spike value is `tau0^2 C(T)`, where `C(T)=sum_h,l a_hl^2`; `E_T C(T)=H_g c_g`. Substituting `E C` inside a reciprocal is a further working approximation. For general leaf-state averaging, sum over the independent states of every occupied leaf. A profile reaches one leaf per tree, but a standardized mean generally reaches many. The two-leaf check in the evidence file gives exact leaf-state average precision 12.3156 versus 20.9615 for the shortcut. Even the main argument that the averaged ESS can never have order n2 is too strong: all-spike probability can be positive and w can approach one.

The reported realized ESS is produced by `forest_gvar_cpp`, which uses the correct squared proportions; `lrc_ess_realized` returns that result as `mean` and returns the c_g approximation separately as `mean_cg_approx`. The wrapper reports `mean`. Therefore the manuscript's approximate realized formula does not describe the results. Nor is V_g at posterior states a posterior variance: the posterior also has leaf covariance and mean/mixture uncertainty. It is a conditional prior-variance diagnostic evaluated at posterior states.

**Concrete fix.** Define a_hl, C(T), and V_g; display the general leaf-state formula (or omit the unused full-mixture formula); state exactly where the implemented all-spike calibration plugs in H_g c_g. Define realized ESS as posterior averaging of `sigma1^2/(V_mu^f+V_g(T,z,tau))`, with V_mu^f fixed from Stage 1. Preserve the regional definitions and state this is a diagnostic, not an identified number of actually borrowed patients. Restrict the Binomial simplification to pointwise/stump cases or an explicitly different tree-wide-state model. Remove its claimed implication about all order-n2 targets.

**Downstream work.** The used calibration and correct realized output can remain unchanged under the above clarification. The exported `lrc_ess_prior` helper still implements the tree-wide-state shortcut and should be documented as such or corrected separately. A repository call-site search found no result-generating call to that helper. Changing c_g plug-in calibration to exact partition averaging would change fitted scales and require refits; do not silently make that change.

## M3. Reconcile Appendix D's theorem, block limit, ceiling anchor, and actual boundary rule

**Location.** Main Theorem 3 and discussion (`sections_1_2.tex:185-194`); Appendix D (`appendix_new.tex:399-443`); `lrcbart.R:403-451`.

**Issue and consequence.** Three claims overstate what is proved or implemented.

1. At fixed block length m, the theorem defines G_m using the raw-block law, but the empirical rule rescales by gamma_B. Its limit is gamma U, with gamma=V_*/E U and U the raw block variance, generally gamma != 1 under autocorrelation. The proof acknowledges this and then incorrectly calls the result G_m of the raw law. The CLT is circularly justified by “Theorem 3(c) of the paper,” while the necessary moments and stochastic equicontinuity for the random-threshold empirical process are not supplied. Geometric ergodicity alone does not provide all required CLTs for unbounded variance statistics.
2. The cap branch tests the whole-chain ceiling and returns the smallest grid value regardless of G(s_min). Thus it is a separate override, not a consequence of the block-quantile Pr-rule or its population positivity criterion. The largest-grid fallback can also return a value failing q. A small spike scale implies complete pooling only in the all-spike model w=1; a robust mixture with slabs remains a robust mixture.
3. The fixed single-tree finite-prior anchor is a lower bound on the precision ceiling; it does not imply that the ceiling exactly equals scaled n2 under matched proportions, or that every covariate shift lowers it. Those statements hold in the flat-prior limit with RWD support in all relevant leaves. With n2=100, sigma2^2=1, lambda_f^2=.1, and pi=(.9,.1), changing a from pi to (5/6,1/6) reduces V from .0086 to .008333 and raises the ceiling from 116.28 to 120. Also H_f lambda_f^2 is constant under the stated range rule, so H_f=50 alone cannot justify negligible prior information.

**Concrete fix.** Keep the elementary monotonicity theorem for the stated working target. For fixed m define the population law of the rescaled statistic; assume finite stationary second moment, positive finite E U, and no probability at the indicator threshold at queried scales. An ergodic law and a sandwich argument give the pointwise result. A finite-grid selection statement also needs a nonempty admissible set and no G_m=q grid point. For fixed B and increasing m, state convergence away from a grid exactly equal to the root. Remove the CLT/rate unless a complete independent proof with explicit assumptions is supplied. Separate mathematical exact-integral statements from finite inner Monte Carlo and c_g estimation. Describe the whole-chain cap and upper-grid failure branches separately; neither branch certifies the Pr-rule inequality. Restrict exact sample-size/overlap equalities to the flat-prior single-tree anchor, require pi_l>0 whenever a_l>0, and explain that a zero-RWD leaf retains prior variance at finite lambda_f. State the positive floor and its feasibility condition when invoking the truncated law. Remove the unsupported H_f-based negligibility assertion.

**Downstream work.** This repair can preserve numerical outputs by describing the implemented override and working target. Altering the cap, truncation law, or finite-grid selection would require recalibration and refits. The theorem is a conditional mathematical result and must not be presented as evidence that the actual BART chains meet its assumptions.

## M4. Distinguish population nonidentification from posterior extrapolation through shared leaves

**Location.** `sections_1_2.tex:93`; Proposition E (`appendix_new.tex:449-465`); single-arm/support discussion in Sections 2.6 and 3.4.

**Issue and consequence.** “Where the RCT has no support g is a prior draw” and Proposition E(d)'s assertion that its law is the propagated prior are false for the tree posterior in general. The finite-sample remark immediately supplies the counterexample: a leaf can straddle supported and unsupported regions, so its RCT-informed parameter predicts both. Learned partitions and global hyperparameters also transmit information. Nonparametric population identification asks what source-specific regression functions are determined without extrapolation restrictions; it is different from whether a fitted finite tree extends a learned value beyond support. There is also no population zero-support region merely because a smooth logistic assignment makes RWD sparse: the Sc2 design has poor finite-sample overlap.

**Concrete fix.** Retain `f` identified on external support, `f+g` on RCT-control support, and their difference on the intersection as unrestricted population statements. Outside the intersection say the unsupported component is not identified from source-specific outcome data; posterior values depend on tree extrapolation and priors. The exact prior conditional claim is leaf-wise: an unobserved g leaf, conditional on partitions and shared hyperparameters, retains its leaf prior. In the true single-arm case the entire g process is independent of the likelihood and remains prior-distributed. State variance formulas conditional on shared scales and w, or average over those parameters. Use “limited RWD support/overlap” in the Sc2 finite-sample interpretation.

**Downstream work.** No refits are needed for these scope corrections. New assertions about full-posterior identification or extrapolation consistency require additional theory rather than relabeling the current proposition.

## M5. The G.6 optimized cost and its posterior conclusions are not established

**Location.** `web_appendix_G.tex:125-127`; `sections_1_2.tex:53,93` and associated design claims.

**Issue and consequence.** The cost `m c + delta^2/(2 m lambda^2)` treats the remaining spike contributions as zero and its quoted threshold comes from unconstrained continuous minimization over m. In an ensemble, m is an integer in 1,...,H, so the continuous optimum may be infeasible. More fundamentally, after fixing which m leaves are slabs, all H contributions can share the shift. Minimizing the Gaussian quadratic penalty subject to their sum being delta gives

`Q_m(delta) = m c + delta^2 / [2{(H-m) tau^2 + m lambda^2}], m=0,...,H`.

This follows from Lagrange multipliers: each contribution is delta times its variance divided by total variance. For c>0 and lambda^2>tau^2, comparison with Q_0 shows that the first possible improvement is at m=1 with threshold squared `2 c H tau^2 {H tau^2 + lambda^2-tau^2}/(lambda^2-tau^2)`. A concrete check (H=5,tau^2=.01,lambda^2=1,c=1,delta=.2) satisfies the manuscript's threshold but has Q_0=.4 and Q_1=1.0192, so no slab allocation improves the joint prior-mode objective.

Even a correct conditional prior-mode comparison does not determine posterior indicators. Integrating leaf means introduces determinant factors and configuration multiplicities; the likelihood and candidate tree partitions matter. The assertion that f's compromise pays in proportion to n1+n2 also omits minimization: a Gaussian common-mean compromise has effective information proportional to A1 A2/(A1+A2), where Ag=ng/sigma_g^2, and is limited by the weaker source. The numerical 17–86 threshold and the absolute statements about clinically sized shifts therefore do not follow.

**Concrete fix.** Retain G.6 only as an explicitly conditional illustration of a prior penalty, with a constrained discrete objective, stated assumptions and a short derivation. Remove the 17–86 values unless their complete inputs and corrected calculation are verified. Delete the general impossibility claim for the shared-leaf design and the assertion that the proposed model necessarily puts a shift in one or two leaves. Say separate ensembles permit different priors and partitions, with the actual performance supported by the simulations. Avoid equating a prior mode with an indicator probability or mixing behavior.

**Downstream work.** Narrowing this motivation needs no numerical replacement. A comparative posterior claim about alternative architectures needs a properly specified comparator, fitted analyses, and diagnostics.

# Minor comments

## m1. Match pointwise summaries to their distinct producing calculations

**Location.** Section 2.4 final sentence, Figure 2/application ESS caption; `lrcbart.R:514-518` and `02-validation/sc4_boundary/pointwise_ess.R:68-73`.

**Issue and consequence.** The application map uses `sigma1^2/[V_f(x)+H_g mean(tau0^2)]`, while the secondary validation helper averages the reciprocal over prior scale draws. These are different by Jensen's inequality. An unqualified “prior pointwise ESS” makes the reported map and supplemental values appear to use the same integration convention.

**Fix.** Define the application map as a plug-in working precision evaluated at the prior mean spike variance (its posterior estimate in single-arm mode), and state the integrated convention for the validation summaries. Preserve the computed numbers and clarify rather than silently recomputing them.

## m2. Add explicit interior-w conditions and retain the comparison basis

**Location.** Main Theorems 1–2 and the subsequent interpretation (`sections_1_2.tex:151-182`); Appendix B.

**Issue and consequence.** The main statements invoke log odds and finite C under Assumption A1 alone, but these quantities need 0<w<1. Also the theorem bounds the difference from a slab-regularized leaf, not an unregularized RCT-only estimate or a frequentist bias under the full fitted model.

**Fix.** Add 0<w<1 to the conditional results using those expressions, handle w=0/1 as separately stated limiting cases, and keep the exact slab comparator in nearby prose. Preserve the existing full-posterior limitations.

## m3. Restrict claims that prior discrepancy noise leaves point estimates unchanged to linear means

**Location.** `sections_3_5.tex:101,127-133`; G.7 `web_appendix_G.tex:137-139`; `tables/tab_single_full.tex`, `tables/tab_app.tex`.

**Issue and consequence.** Zero-mean g implies no change in the expected control log-time mean, but not in nonlinear survival summaries. The full single-arm table already shows median-ratio bias .076 for the pooled row versus .108 at target 100 under Sc1,nT=200, and much larger changes for w=.9. The application table gives PFS mean .97 at w=1 versus 1.00 at w=.9. Thus “every design ... within .005” and “the mean unchanged” exceed the preserved numbers.

**Fix.** Limit exact mean-preservation statements to the Gaussian/control log-time linear mean. Describe nonlinear survival means as empirically changed by the prior draw, retaining the existing table entries. When “nearly unchanged” is used, identify the actual target ladder/estimand over which it is supported.

## m4. Keep exact-kernel invariance separate from finite-chain independence

**Location.** Appendix A proof at `appendix_new.tex:107` (“then yields a draw from the joint conditional”).

**Issue and consequence.** One MH move followed by a conditional leaf draw defines a kernel preserving the joint conditional; it does not yield an independent draw from that distribution from arbitrary initialization.

**Fix.** Say the composite update leaves the joint conditional invariant. Preserve the proposition's mathematical conclusion and avoid treating it as a mixing or implementation-certification claim.

## m5. Remove internal production references from scientific proofs when revising them

**Location.** `appendix_new.tex:110-111,312-316,377-378,427,431,436-437,468-470`.

**Issue and consequence.** References to “the specification,” an unspecified “verification report,” and superseded “Theorem 3(c) of the paper” make the supplement depend on undocumented production history; the last is also mathematically circular as addressed in M3.

**Fix.** Give a self-contained statement of current assumptions, derivations and limitations. Keep historical revision accounting in the review ledger. Do not cut substantive caveats to achieve brevity.

# Suggestions to enhance the paper

A short reproducibility note distinguishing (i) the fixed-partition conditional variance, (ii) the all-spike calibration approximation, (iii) the realized posterior-state diagnostic, and (iv) the pointwise plug-in map would resolve much confusion without a longer theoretical agenda. A full-model empirical reliability study would be useful future work, but is not a prerequisite for honestly reporting the current simulations.

# Source-dependent flags

- This review did not reproduce patient-level preprocessing or reassess clinical validity of conditioning on post-induction ASCT. The current text itself labels that variable as baseline; the contextual reviewer was asked to assess the application's estimand and interpretation. A claimed causal treatment benefit would need separate clinical/causal qualification.
- The earlier MAP-AFT-BART sampler is stated to be unavailable. Attribution of the earlier extreme tails to a particular implementation mechanism cannot be independently verified from the current results alone; retain a model-level explanation only, not a proven decomposition of numerical causes.
- No full simulation or application fit was rerun. The algebra checks concern definitions and conditional objectives, not operating characteristics or sampler convergence. All reported numerical results remain the existing outputs.

# Verdict

The core modeling construction is reviewable and the conditional robustness results are useful. The necessary revision is to make the calibration target, leaf aggregation, theorem conditions, implementation exceptions, extrapolation, and architecture motivation scientifically consistent. The existing results need not be silently replaced to do this. If an exact-ELIR or different calibration method is adopted, a new empirical round is required before its results can replace the preserved evidence.

VERDICT: major revision
MAJOR COUNT: 5

What bothered me most: the paper currently treats several different information summaries as one exact ELIR while its correct leaf-based realized calculation is more careful than its prose.


---

# Summary

Independent applied/contextual review of the complete current manuscript assembly, both main-text files, Appendices A–F and G, the complete `main_lrcbart.tex` deck, and its model-diagram insert. I did not read the methodologist's review or edit manuscript, code, figures, or slides. Source checks below use the producing code/CSV where available and the primary CAHB paper. This is an independent AI review, not external human approval.

The paper has a useful, concrete contribution: a small discrepancy ensemble separates the regularization of source differences from that of the external-control outcome surface, and the output distinguishes compatibility, control-surface accuracy, and information summaries. The current text nevertheless overstates the shared-leaf comparison, misdescribes a closest comparator as parametric, confuses population nonidentification with an unchanged posterior prior, and carries incorrect single-arm claims onto nonlinear survival scales. The application also needs an explicit descriptive interpretation because ASCT after induction is included as though it were a baseline covariate. These issues can largely be corrected without changing any reported numerical result.

# What the paper does well

- The external-centered explanation in Section 2.2 and slides 5–8 is clear and coherent. The diagram correctly distinguishes one leaf contribution from the whole discrepancy, and identifies the shared scales and weight.
- Section 3 distinguishes all-profile surface error, control-profile checks, regional compatibility, and regional realized ESS. The far-boundary analysis prevents the surface loss from being dismissed as only a boundary effect.
- The manuscript now describes CAHB as a transfer with explicit smoothing, covariate-selection, precision-scaling, and joint-draw choices and reports its variants. Its own precision share is appropriately separated from LRC's compatibility probability.
- The fixed-partition and fixed-surface limits of the map results, the inability of a single-arm trial to learn control discrepancies, and the lack of Type I error control from prior-ESS calibration are stated in useful places.
- The result tables and comparison code make several remaining prose errors directly checkable without rerunning simulations.

# Major comments

## M1. Remove the unsupported impossibility and indicator claims from the shared-leaf comparison

1. **Location:** `sections_1_2.tex:44,93`; G.6 (`web_appendix_G.tex:127`); deck “Why a separate discrepancy ensemble” (`main_lrcbart.tex:218–224`). Quotes include “which a single ensemble with two means per leaf cannot do,” “no shift of clinical size flips an indicator,” and “the spike cannot hide it by spreading it.”
2. **Issue and consequence:** The calculation compares selected optimized prior penalties, then treats their ordering as a result about fitted posterior indicator probabilities and detection. It omits partition/proposal possibilities, integrated probability mass, other allocations of the shift across spike and slab leaves, and likelihood information. Optimizing `m c + delta^2/(2 m lambda^2)` continuously also requires `c>0` and an optimum within the allowed integer range; the displayed threshold is not an unconditional switching threshold. A shared partition can in principle represent regional differences, so the claimed impossibility overstates the methodological distinction. “Clinical size” is not defined by a simulated outcome range.
3. **Concrete fix:** Keep G.6 only as a carefully stated illustrative prior-penalty comparison with explicit fixed partitions/scales and no inference about posterior switching. The main text and slide should say that separate ensembles permit different partition and shrinkage priors for the baseline surface and discrepancy; keeping `H_g` small limits ways of distributing a discrepancy under the spike. Remove the unsupported 17–86 clinical threshold and guaranteed flat-map claims unless a fully specified derivation and relevant empirical comparison actually establish them. Do not replace them with another unverified superiority claim. No simulation results need change.

## M2. Correct CAHB's nonparametric classification and narrow the novelty claim

1. **Location:** `sections_1_2.tex:49,51`, especially “closest parametric relative” and “neither on a nonparametric outcome surface”; G.0 (`web_appendix_G.tex:15`) repeats the former.
2. **Issue and consequence:** CAHB explicitly estimates the current means and local precision nonparametrically. The verified primary paper states this in Section 2 (journal p. 5341, PDF p. 4) and Section 4 (journal p. 5343, PDF p. 6). Its historical model in one example does not make the method's outcome estimation parametric. The current exclusion materially inflates the novelty claim and undermines the fairness of the closest-comparator discussion. Source: [Jin et al. (2023), published paper](https://jinhuaqing.github.io/files/2023_SIM_CAHB.pdf), DOI 10.1002/sim.9913.
3. **Concrete fix:** Call CAHB the closest kernel-local borrowing relative. State that it estimates covariate-dependent current mean and precision functions through local kernel regression. Distinguish LRC's learned tree partitions, jointly estimated external surface and discrepancy, and explicit leaf mixture/calibration, without claiming that earlier local borrowing lacks a nonparametric outcome surface. Keep the existing transfer caveats, paired `R=200`, and variant results. Use the paper's term “covariate-adaptive historical control borrowing” consistently instead of expanding CAHB as “hybrid borrowing.”

## M3. Separate nonidentification from prior or extrapolated prediction off support

1. **Location:** `sections_1_2.tex:93`; Appendix E, Proposition “Population identification,” part (d), `appendix_new.tex:455`; compare the finite-sample qualification at line 465.
2. **Issue and consequence:** Lack of RCT support at a profile establishes that the discrepancy is not identified there from a local outcome distribution. It does not make its full posterior an independent prior draw. A tree leaf can contain RCT controls elsewhere and extend their fitted discrepancy to that profile; partitions and shared hyperparameters are also learned. The appendix's finite-sample remark already gives this counterexample, so the stronger proposition and main-text sentence contradict the model's own description. On common support, separation of response surfaces follows from the two conditional means in the population; shrinkage helps finite-sample regularization and is not itself identification.
3. **Concrete fix:** State the population identification conclusions only. Say that outside common support, separation relies on the model's extrapolation and prior assumptions. A reached leaf with no RCT observations has its prior conditional at fixed partitions and current hyperparameters; a leaf spanning observed and unobserved support extrapolates a fitted discrepancy. Reserve “g remains a prior draw everywhere” for the true single-arm model with no RCT controls anywhere. Keep the main-text correction short and reconcile the proposition/proof with the existing finite-sample remark.

## M4. The application needs a descriptive estimand and correct ASCT timing

1. **Location:** `sections_3_5.tex:113` calls all six variables baseline covariates while including “receipt of ... ASCT after induction”; lines 120–133 discuss standardized treatment benefit. G.5 table caption at line 80 and deck “EloKRd and UCMM data,” line 434, also label them baseline.
2. **Issue and consequence:** The implemented `asct` variable is eventual recorded transplant receipt after time zero, not a baseline attribute. `04-application/code/prep_data.R:15,33–41` includes `transplant_off_protocol`; original UCMM preparation excludes transplants on/before induction and then defines the flag from whether a transplant date exists (`../yunxuan-repo/mapbart-case-study-mm/data_cleaning_ucmm.R:61–62,159–160,207`). Standardizing predictions at treated patients' post-induction ASCT statuses does not by itself identify a total causal treatment effect. The single-arm design also cannot learn the control discrepancy, even when registry support is ample. Calling the analysis evidence for/against marginal treatment benefit without this qualification invites an unsupported clinical interpretation.
3. **Concrete fix:** Preserve the actual fitted covariates and all numbers. Replace “six baseline covariates” with “five baseline covariates and recorded ASCT after induction,” and use “cohort characteristics” for the table caption. Define the application as a model-based standardized comparison conditional on these recorded profiles and the chosen discrepancy prior. State briefly that ASCT occurs after induction, the adjustment does not establish a total causal effect, and no internal control can verify residual exchangeability. Revise conclusions to describe the fitted RMST comparisons, without an efficacy or inferiority conclusion. A baseline-only or time-aware causal analysis is separate author/analysis work; do not imply it was performed.

## M5. Zero-mean discrepancy does not preserve nonlinear survival-ratio means

1. **Location:** `sections_3_5.tex:101,131`; G.7 (`web_appendix_G.tex:137,139`); deck “The single-arm design,” line 413. Quotes: “In every design ... within 0.005,” “with the mean unchanged,” and “point estimate unchanged.”
2. **Issue and consequence:** Zero-mean `g` preserves the expected standardized Gaussian response mean under the stated independence, but survival quantiles, RMSTs, and ratios are nonlinear. The code transforms each control draw before taking the ratio (`01-code/lrcbart/R/lrcbart.R:554–570`), and `metrics.R:53` summarizes by its mean. The stored results contradict the claimed invariance. In `02-validation/full/single_arm/summary_table.csv`, Sc1, survival `n_T=200`, median-ratio bias is 0.108027 for target-100 `w=1`, 0.076269 for `LRC-BART-full`, and 1.491964 for target-100 `w=0.9`. At `n_T=30`, the `w=1` and `w=0.9` biases are 0.180291 and 1.423366. In `04-application/results/lrcbart_summary.csv`, harmonized PFS RMST means at target 30 are 0.973861 (`w=1`) and 0.997027 (`w=0.9`); the medians are much closer, 0.964111 and 0.966127.
3. **Concrete fix:** Restrict exact mean-preservation reasoning to the latent response/log-time mean. Describe near-preservation of Gaussian point estimates as empirical, and explicitly state that survival-ratio means may increase substantially when discrepancy variance increases. Replace the application's “mean unchanged” with the verified small observed change, or simply say the interval widens and the estimate is similar. Correct G.7 and the slide accordingly, retaining every existing table result. Do not substitute posterior medians while claiming the same estimand summary.

# Minor comments

## m1. Synchronize the deck's quantitative interpretation and qualifications

1. **Location:** Full deck, especially lines 246–277, 360–368, 475–501; exact audit below.
2. **Issue and consequence:** Several statements predate corrections already present in the paper: one residual SD instead of two thirds; an established sample-size limit; additive ceilings; no qualification on the single-tree flat-prior sample-size anchor; and a blanket interpretation of the compatibility map as information used. The reader would receive a stronger and in places false account from the talk.
3. **Concrete fix:** Apply the scoped deck audit below without changing the 25-slide structure or inventing numerical evidence. Reuse the manuscript's verified regional ESS and far-boundary summaries in existing blocks. Keep the centered-prior slides and model diagram intact.

## m2. State the realized-ESS and calibration units as working information summaries

1. **Location:** `sections_1_2.tex:122–141`; Appendix D; application at `sections_3_5.tex:122,133`.
2. **Issue and consequence:** The application correctly calibrates the control log-time mean, but the rest of the prose repeatedly reads as if a variance-equivalent scalar ESS were exact marginal ELIR or a count of external patients used. A standardized ensemble can span many leaves, so its variance depends on leaf proportions and states, not a pointwise count of reached leaves. The local ESS is not additive across profiles or regions and does not itself demonstrate compatibility.
3. **Concrete fix:** Once the producing calibration audit is resolved, label the implemented criterion consistently as a working inverse-variance/conditional-information summary, with its approximation and conditioning explicit. Retain the statement that application calibration concerns log-time mean rather than RMST, and retain the regional nonadditivity caveat. Avoid replacing a code-dependent claim with a new formula without tracing its implementation and numerical implications. This comment is scoped to interpretation; the parent is auditing the equations and calibration implementation.

## m3. Correct the coding mechanism without claiming identical predictions for every model

1. **Location:** G.5 (`web_appendix_G.tex:75,109`) and deck lines 439–440: “every control model ... predicts ... from the 13-patient Hispanic stratum.”
2. **Issue and consequence:** `prep_covariates_orig` creates two different non-Hispanic dummy columns. The trial's retained registry non-Hispanic column is zero, matching Hispanic coding on that variable, but a trial-only dummy has no observed variation among controls. In the Stan reference fit the corresponding coefficient is retained with prior uncertainty (`run_reference.R:77–85,88–100`); predictions also depend on other covariates and model structure. The identified coding defect is real, but the current universal prediction statement is too literal.
3. **Concrete fix:** Explain the two incompatible dummy columns and the unsupported trial profile coding. Say the registry non-Hispanic indicator takes the Hispanic-coded value for trial patients, with model-dependent extrapolation and prior effects. Keep all coding-sensitivity numbers; do not claim the coding explains a fixed fraction of every model difference.

## m4. Align cohort restriction prose with the producing script

1. **Location:** `sections_3_5.tex:113` says restrictions occur “in order” beginning with recorded regimen then positive survival; G.5 line 105 uses another order.
2. **Issue and consequence:** The UCMM script first excludes negative OS/PFS times, then requires a nonmissing induction start, then excludes ASCT on/before induction, then E-Rd/Elo-Rd (`data_cleaning_ucmm.R:58–77`). It does not implement a generic “positive survival times” or “recorded induction regimen” rule at those steps. A reader cannot reproduce the cohort from the current description exactly.
3. **Concrete fix:** Either describe the four actual masks in their implemented order or omit “in order” and give a concise accurate cohort definition. Preserve 797/253 and the existing analytic cohort. No cohort rebuild is required for this textual correction.

## m5. Remove residual production history and repeated narration to meet the length target

1. **Location:** Main text `sections_1_2.tex:38,42–53,89,98–115,153–195`; application `sections_3_5.tex:129–133`; G.5 lines 75/109 repeat the coding account. Appendix A, C, D and F contain “previous model,” “specification,” “verification report,” and “Gap” narration.
2. **Issue and consequence:** The main-text target is exceeded and several paragraphs repeat material already displayed mathematically or placed in the appendix. The historical/production explanations make the paper read as a revision record. Some Appendix D references to old Theorem 3(c) are not self-contained current proofs.
3. **Concrete fix:** Cut approximately 600–800 words from the main text after the scientific corrections: shorten the study-level literature tutorial (~120), merge contribution/roadmap repeats (~90), move conjugate-update detail already in Appendix A out of prose (~130), compress theorem retellings while retaining assumptions (~140), and shorten repeated coding/earlier-tail application history (~120). Treat these as approximate targets and recount. Retain CAHB choices, support restrictions, fixed-partition conditions, nonlinear-ratio limitations, and empirical surface losses. Since Appendix A–F is protected, make only justified scientific/cross-reference corrections there and log any wider cleanup as author work rather than silently rewriting it.

## m6. Tighten small literal errors without changing substance

1. **Location:** `sections_1_2.tex:89` has a comma before “The sampler”; line 117 uses `log y_i` although line 63 already defines survival `y_i` as log time; `sections_3_5.tex:49` says C-BART “cannot let the sources differ”; application diagnostic statement line 120 says RMST ratio whereas code line 75 assesses `log(m)`.
2. **Issue and consequence:** The punctuation and outcome notation interrupt reading; C-BART permits a discrepancy even without a slab; the scale used for diagnostics is reproducibility information.
3. **Concrete fix:** Correct the comma, use a distinct raw survival-time symbol consistently, describe C-BART as having one global commensurability scale without local slab adaptation, and specify that split-Rhat/bulk ESS were assessed on the log RMST ratio. Preserve diagnostic values.

# Suggestions to enhance the paper

The following is the requested exact stale-deck audit, supporting M1–M5 and m1–m3. Locations are current `main_lrcbart.tex` line numbers; titles remain unchanged.

| Location / current claim | Scoped correction |
|---|---|
| Subtitle, line 64: “Borrow where ... abstain where ...” | Frame as an aim (“Model local compatibility and calibrate borrowing”) or keep explicitly aspirational; avoid a guarantee contradicted by moderate shifts. |
| “Why a separate discrepancy ensemble,” 218–224: thresholds 17–86, clinically sized shifts cannot flip indicators | Replace the block with the qualified distinct-regularization rationale in M1; no new comparative performance claim. |
| “Every leaf marginal is exact,” 239: `tau_0^2` and `w` update “from the spike leaves” | `tau_0^2` uses spike-leaf means/count; `w` uses spike and total leaf counts. Say “from the leaf states and means.” |
| “Prior effective sample size,” 248: calls the variance expression expected ELIR | Match the verified working-ESS definition and its Gaussian/conditioning limits. |
| Same slide, 256: matched distributions imply sample-size equality; shift necessarily lowers it | State the single fixed-tree, flat-prior, occupied-leaf anchor. The fitted BART ceiling is computed from its posterior variance. |
| “Three results,” 274: map limit without representability and threshold-boundary exception | Add fixed partitions, representable discrepancy, fixed `f`, and `|d| != c`; can be one short line in the theorem block. |
| Same slide, 277: convexity/zero limit without untruncated qualifier; positive solution exactly below ceiling | Specify the untruncated working calibration, finite-grid feasibility, block criterion and Monte Carlo conditions. Mention positive floor if describing the truncated sampling law. |
| “The borrowing map,” 290: all-spike probability “does not settle” | Say it need not identify compatibility as sample size grows at a fixed positive spike scale; it can converge to an interior probability. |
| “Simulation design,” 305/322: only `R=500` and old comparator list | Add CAHB transfer on the first 200 paired Gaussian replicates, with transfer choices/variants described in supplement. Remove “new” from Sc4/Sc5. No new slide or mixed-denominator table required. |
| “What the Gaussian study shows,” 362: “and more on the control surface” | State that surface comparisons depend on evaluation profiles, or delete this extra clause. Existing control-profile Sc1 RMSE is 0.704824 vs 0.699068 (`sc4_boundary/paired.csv`), so it is not a universal surface gain. |
| Same slide, 365: bias is “the price” of the spike | Say the retained spike is a possible explanation; the leaf theorem does not isolate the cause of full-ensemble bias. |
| Same slide, 368: “a sample-size limit that tuning did not remove” | “Underestimation persisted in the examined configurations”; do not infer a demonstrated sample-size limit. |
| Same slide, 368: “At one residual SD, near d_half” for shift 1 | “At two thirds of a residual SD” (1/1.5); the manuscript's d_half is about 0.5, so shift 1 is about twice that illustrative leaf threshold. |
| Same slide, 368: Sc4 surface figures unlabeled by shift; boundary treated as explanation | Tie 0.824/0.785 to shift 2 and all profiles in the compatible region. Include the existing far-control-profile loss, 0.747/0.723, SE 0.005, `R=200`, by shortening surrounding narration. |
| Same slide or “The borrowing map”: regional realized ESS omitted | Incorporate the existing shift-2 regional values 9.5/0.44 (`R=200`) and say the two regional estimands are not additive. A compatibility probability is not this ESS. |
| “Survival outcomes,” 393: “borrowing makes it usable” and “through the convex ratio” | State the observed RMSE reduction and that nonlinear transformation and wider posterior contribute to mean bias; a general ratio is not jointly convex. |
| “The single-arm design,” 413: slab leaves point estimate unchanged | Limit exact mean-centering to latent log-time; survival ratio means can change substantially (M5). |
| “EloKRd and UCMM data,” 434/440 | Correct post-induction ASCT timing and explain incompatible dummy coding without asserting identical Hispanic-stratum predictions for every model (M4, m3). |
| “EloKRd results,” 467: “No ... supports marginal benefit”; prior proven cause of earlier tails | Use the conditional descriptive comparison in M4. Earlier sampler unavailable: report observed shorter tails without a verified attribution to that prior alone. |
| “The local ESS map,” 476: target 30 “spread thinly” | Say target 30 concerns the standardized log-time mean and does not imply equal information at each profile. Local values are not allocations whose sum is 30. |
| “Summary,” 483/486–487: tree prior chooses borrowing; unconditional guarantees | Say posterior leaf states vary over learned partitions; compatibility map reports discrepancy probability. Retain conditional theorem limits and Sc4 surface/coverage losses. |
| “Limitations,” 493: posterior over partitions is “where the boundary leak lives” | Remove the unsupported mechanistic attribution; state the absence of full-posterior theory and observed surface loss. |
| Same slide, 495: several sources, “the ceiling becomes a sum” | “Source-specific discrepancies would need joint calibration; additive ceilings have not been established.” |
| Takeaway, 501: for every profile it reports how much external information “it used” | “Reports local compatibility and information summaries”; explicitly distinguish the prior pointwise ESS from realized regional ESS and the single-arm absence of learned compatibility. |

The model-diagram insert (`inserts/lrc_model_diagram.tex`) needs no scientific correction on the issues reviewed: its external centering, per-profile sum, and shared scales/weight are consistent with the model.

# Source-dependent flags

- **Earlier MAP-AFT-BART tail attribution:** G.5 line 109 and deck line 467 assert the previous prior caused the earlier long tails, but `run_reference.R:6–9` states that the former sampler is unavailable and these rows were not rerun. The new distributions and old printed intervals are observable; a controlled attribution among prior, model, coding and implementation is not verified. Narrow to observed comparisons; do not call the cause established.
- **Author field:** `sections_1_2.tex:19` is empty. Presenter and collaborator credits are not enough to infer final author order, affiliations or corresponding-author metadata. Preserve the field pending author-supplied information; do not fabricate it.
- **Clinical provenance:** The dataset/code verifies the recorded cohort and analysis variables; it does not independently establish the clinical definition, timing, or scientific suitability of every endpoint/covariate or the complete trial/publication metadata. Any new clinical claim or clinical reference needs the actual protocol/publication or data dictionary. No new clinical analysis was run for this review.
- **Calibration equations:** Exact marginal ELIR and multi-leaf aggregation require the separate mathematical/code audit. This review does not certify them or authorize relabeling an approximate criterion as exact. All recommended interpretation changes must follow that audit and preserve numerical lineage.

# Verdict

The central construction and simulation record remain useful. Five blocking interpretation/positioning errors need correction, with a synchronized deck and a shorter, self-contained main text. The remedies above preserve the results and do not require inventing new analyses. The application should remain a descriptive illustration unless a separately specified causal analysis is performed.

VERDICT: major revision
MAJOR COUNT: 5

What bothered me most: several carefully qualified result paragraphs are contradicted by broader methodological and single-arm claims elsewhere in the paper and talk.
