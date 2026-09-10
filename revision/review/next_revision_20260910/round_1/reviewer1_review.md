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
