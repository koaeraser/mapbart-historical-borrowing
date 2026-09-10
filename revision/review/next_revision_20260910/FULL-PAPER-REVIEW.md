# Full-paper review: LRC-BART

Review date: September 10, 2026. Saved at the user's request from the side conversation.

**Overall assessment: major revision.** The Gaussian leaf integrals, mixture conditionals, conditional influence bound, and fixed-partition map-consistency arguments appear coherent. The principal problems concern their interpretation, the calibration framework, several claims contradicted by the numerical tables, and the application.

## Scope and reviewed version

This was a read-only review by the side-conversation assistant, separate from the main thread's revision work. It is an internal assessment, not external peer review or scientific approval.

The review covered Sections 1–5, Web Appendices A–G, relevant calibration and sampler code, and selected primary references. All 82 pages of the clean PDF were inspected at overview scale, with key figures and tables inspected more closely. No manuscript, code, or results were changed, and simulations were not rerun. This file is the only artifact created by the side conversation.

The reviewed PDF was 05-writing/build/main_draft.pdf, with a creation time of September 10, 2026, 09:28:57 CDT. The sources read were:

- 05-writing/sections_1_2.tex
- 05-writing/sections_3_5.tex
- 03-theory/appendix_new.tex
- 05-writing/web_appendix_G.tex
- Main and supplementary tables in 05-writing/tables/

Locations below refer to the version read during this review. The main thread may subsequently have changed these files. Reconcile the findings against the current candidate before treating them as still open.

## Major findings

### 1. The ELIR identification is incorrect

**Location:** Section 2.4, especially sections_1_2.tex lines 122–141; Appendix D, including the units-and-ceiling remark.

Section 2.4 equates expected local-information-ratio ESS with sigma_1^2 / V for an arbitrary scalar prior with variance V. This identity holds for a Gaussian prior under the specified Gaussian information unit, but generally fails for mixtures and other nonnormal priors. Averaging component information also differs from computing information in the marginal mixture.

The cited Neuenschwander et al. paper explicitly distinguishes these quantities. Under constant unit information 1 / sigma_1^2, marginal ELIR uses the prior expectation of the negative second derivative of the log marginal prior density, multiplied by sigma_1^2. A variance alone does not determine this curvature.

A numerical illustration checked during this review used the equally weighted mixture N(0,1) and N(0,9), with unit residual variance. Its inverse-variance quantity is 0.2000, its marginal information is approximately 0.3455, and its average component information is 0.5556. These are three different quantities.

**Required correction:** Decide whether to retain the current calculation as an explicitly defined variance-based calibration quantity or adopt marginal ELIR. The latter would require recalibration and assessment of downstream results. Renaming alone must not retain unsupported claims about exact marginal information or predictive consistency.

Primary source: [Neuenschwander et al., Predictively Consistent Prior Effective Sample Sizes](https://www.tonyohagan.co.uk/academic/pdf/ess.pdf), Sections 2.4–2.7.

### 2. The standardized-mean ESS formula does not respect independent leaf indicators

**Location:** Section 2.4, sections_1_2.tex lines 122–141; Appendix D line 387; 01-code/lrcbart/R/lrcbart.R lines 466–508.

For fixed partitions, let a_(h,l) be the proportion of target profiles in leaf (h,l). The conditional prior variance of the standardized discrepancy is

\[
V_g=\sum_{h,\ell}a_{h\ell}^{\,2}
\{z_{h\ell}\tau_0^2+(1-z_{h\ell})\tau_1^2\}.
\]

The manuscript's binomial count over H_g trees does not generally represent this distribution when the standardized estimand spans multiple independently labelled leaves. A pointwise profile reaches one leaf per tree; a standardized mean can span many leaves per tree.

For example, with one tree, two equally weighted leaves, w=0.5, tau_0^2=1, tau_1^2=9, V_mu^f=1 and sigma_1^2=1, exact averaging over the four leaf-state configurations gives approximately 0.3550. The tree-level binomial expression with c_g=0.5 gives approximately 0.4242.

The realized-ESS implementation already computes the leaf-weighted variance through forest_gvar_cpp and returns its corresponding summary as the primary result. It reports a c_g approximation separately. The manuscript instead prints the approximation as the realized-ESS definition. The pointwise map also uses a plug-in mean scale in the code, which should be distinguished from averaging reciprocal variances.

**Required correction:** Reconcile prior calibration, prior-averaged reporting, realized regional summaries, and pointwise summaries. Specify what is conditioned on and what is averaged over at each step. Distinguish the variance under the conditional prior component from the posterior variance of g after seeing trial outcomes. Quantify any approximation whose scientific interpretation is retained.

### 3. The calibration theorem and implemented selection rule remain inconsistent

**Location:** Main Theorem 3, sections_1_2.tex lines 185–195; Appendix D, appendix_new.tex lines 399–431; lrcbart.R lines 403–451.

There are three distinct problems:

- Appendix D defines the fixed-block population law using raw block variances, although the estimator uses rescaled variances. The scaling need not converge to one. For a stationary AR(1) process with correlation 0.8 and blocks of two draws, the mean sample variance within a block is 0.2 times the stationary variance, so the mean-matching scaling converges to 5. The limiting population function must correspond to the scaled block statistic.
- The main theorem identifies a finite-grid selection with the continuously defined optimum. In general it can equal only the corresponding grid selection. The appendix also invokes a CLT through references to the earlier paper's theorem; the required moment, regularity and rescaling conditions should be stated directly and verified for the claimed result.
- The code's capped-target branch bypasses the stated q=0.95 Pr-rule. With block variances 0.5 and 1.5, whole-chain variance 1, unit residual variance and target 1, only half the blocks satisfy the zero-scale criterion, yet the cap selects the smallest scale. Thus the implemented cap is an additional rule, not a consequence of the displayed Pr-rule.

**Required correction:** State one coherent population target, its grid approximation, its Monte Carlo conditions, and any explicit capping exception. Check the mathematics and the implemented branch separately. A simulation rerun is needed only if the actual selection algorithm is changed; a revised description must accurately identify the algorithm that produced the existing results.

### 4. The ceiling interpretation overstates the theorem

**Location:** Section 2.4, sections_1_2.tex line 133; interpretation of Theorem 3 at line 195; Appendix D, part (d).

Equality of the ceiling with the external sample size in RCT units requires the stated flat-prior limit. With finite regularization, prior information contributes to the ceiling. The displayed inequality does not establish that arbitrary covariate shift lowers the regularized ceiling.

A checked two-leaf counterexample used n_2=100, sigma_1^2=sigma_2^2=1, lambda_f^2=0.01 and external leaf proportions (0.9,0.1). The ceiling is approximately 229.67 at matched target proportions. Changing target proportions to (0.6333,0.3667) increases it to 300. This does not contradict the flat-prior result; it contradicts extending that interpretation to every finite-prior setting.

The assertion that a fractional target is always feasible also needs the truncation floor and finite search grid addressed.

**Required correction:** Keep the flat-prior anchor separate from the regularized BART quantity. State feasible ranges for the calibration law actually used and report unattained finite-grid targets accurately.

### 5. G.6 does not establish its indicator-switching claims

**Location:** Introduction and Section 2.2; Web Appendix G.6, web_appendix_G.tex lines 125–127.

A configuration containing m slab leaves still permits the remaining spike leaves to contribute to the shift. Its optimized Gaussian penalty therefore involves the total variance

\[
(H-m)\tau^2+m\lambda^2,
\]

rather than only m lambda^2. Under a fixed configuration, a corresponding quadratic cost is m c + delta^2 / [2{(H-m)tau^2 + m lambda^2}]. Optimizing a route that forces all remaining spike contributions to zero does not establish the optimum of the full mixed configuration.

Configuration multiplicity, integer constraints and likelihood contributions also matter. A prior-density comparison alone cannot establish posterior indicator probabilities. The conclusion that no clinically sized shift flips an indicator, and the claim that a shared-leaf ensemble cannot support local borrowing, are unsupported by the supplied calculation.

**Required correction:** Re-derive the constrained comparison with explicit scope or narrow it to a clearly labelled structural motivation. Remove the unsupported clinical-threshold and posterior-indicator conclusions. Reconcile the related slide claim in the main revision.

### 6. Lack of local trial support does not imply a prior draw

**Location:** Section 2.2, final paragraph; Appendix E, especially appendix_new.tex lines 449–465.

A discrepancy-tree leaf can extend across supported and unsupported regions, transmitting posterior information to a profile without local trial observations. Learned partitions and hyperparameters transmit information too. Appendix E's population-identification statement and its finite-sample remark contradict one another on this point.

For example, a stump has a common discrepancy parameter at every profile. Trial observations in one part of its domain update that parameter everywhere, including profiles outside trial support.

**Required correction:** Separate nonparametric population identification, posterior extrapolation through shared structure, and the conditional prior law of a genuinely empty leaf. The global single-arm statement has a different justification because no trial-control likelihood contribution exists anywhere.

### 7. Single-arm survival claims contradict both mathematics and tables

**Location:** Section 3.4; Section 4.3, sections_3_5.tex line 131; Web Appendix G.7, especially line 139; Tables 4, 5 and S6.

A zero-mean discrepancy preserves a linear mean, but generally changes the posterior mean of a survival ratio. For example, if g is N(0,V_g), then

\[
E[e^{a-g}]=e^{a+V_g/2}.
\]

The manuscript states that the slab sensitivity leaves the point estimate unchanged. The tables contradict this:

- Web Table S6, Sc1, median ratio, n_T=30: bias is 0.180 for LRC-BART (100) and 1.423 for LRC-BART (w=0.9).
- Table 5, application PFS: the posterior mean is 0.97 at the primary target and 1.00 under w=0.9.

**Required correction:** Restrict mean-preservation arguments to the linear estimands for which they hold. Reconcile every survival statement with the table and distinguish posterior means, posterior medians, and uncertainty. Treat the observed near-stability of particular RMST estimates as an empirical result, not a general consequence of zero-centered g.

### 8. The application needs an explicit causal interpretation

**Location:** Section 4.1, sections_3_5.tex line 113; Web Table S5; 04-application/code/prep_data.R lines 15 and 41.

ASCT after induction is included among baseline adjustment covariates, despite follow-up beginning at treatment or induction. The model matrix includes the transplant variable. Conditioning on a subsequent treatment can introduce selection or mediator-adjustment problems. Actual timing, the intended treatment strategy and the estimand need assessment before interpreting this as a treatment-benefit comparison.

This is a design concern inferred from the manuscript's timing description, not a finding that the direction or magnitude of bias has been established in the patient data. The equal-residual-variance assumption in the single-arm survival prediction should also be explicit, rather than appearing only as a calibration convention.

**Required correction:** Establish which covariates were known at time zero and justify the role of ASCT. Specify the causal question and assumptions, or limit the analysis to an appropriately described model-based association. Any revised adjustment set would require a documented application reanalysis.

Relevant primary methodological source: [Hernan and colleagues on immortal time, selection and target-trial alignment](https://research-information.bris.ac.uk/ws/files/404569243/hernan_immortaltime_11jun24.pdf).

### 9. The novelty comparison mischaracterizes CAHB

**Location:** Related work, sections_1_2.tex lines 49–51; Web Appendix G.0.

Calling CAHB a parametric relative and claiming it lacks a nonparametric outcome surface conflicts with Jin et al.'s primary paper. That paper uses nonparametric kernel estimation for the current-arm mean functions and the covariate-dependent precision.

**Required correction:** Distinguish LRC-BART through the tree construction, joint model, discrepancy prior and calibration. Do not base the distinction on an incorrect parametric/nonparametric contrast. Preserve the caveats about the transferred comparator implementation and its evaluation on 200 paired replicates.

Primary source: [Jin et al., Bayesian adaptive design for covariate-adaptive historical control information borrowing](https://escholarship.org/content/qt99w3q0wm/qt99w3q0wm.pdf), especially Sections 2–4.

## Writing

The main text contains 7,814 text words including the abstract, as counted by texcount on sections_1_2.tex and sections_3_5.tex. This count excludes the separate table inputs and is distinct from counts that add mathematical expressions, headings and captions.

- The introduction repeats the contribution across several paragraphs. State the clinical problem, existing limitation and precise contribution once, then develop them.
- Results prose frequently reproduces table entries. Retain the key comparisons and their uncertainty, with details left in the tables.
- Appendix F and several Gap paragraphs retain revision-history or development narration. Rewrite the retained scientific content as a self-contained exposition after the technical issues are resolved.
- Terminology drifts between compatibility, borrowing and information. Maintain explicit distinctions among posterior compatibility probability, all-spike probability, calibration ESS, realized variance-based summaries and pointwise prior information.
- The author field is empty.
- The survival notation is inconsistent: y_i is defined as log time, but Section 2.3 subsequently writes log y_i. Define time and log time separately.
- Fix the comma before The sampler in the leaf-prior paragraph.
- Avoid language such as unbiased when the evidence is a small empirical bias in a finite simulation. Keep conditional theorem statements distinct from claims about the full fitted model.

## Presentation

The reviewed PDF is readable, and no obvious clipping was found in the inspected pages. Its existing build log reports no warnings or overfull boxes. This review did not rebuild the manuscript.

- Main simulation tables are dense and use small text. More selective main tables with larger type would improve reading; retain full numerical results in supplementary tables.
- Table 2's regional ESS panel is separated onto the following page, weakening its connection to the compatibility and surface summaries.
- Visible hyperlink borders create unnecessary red and green visual clutter.
- The flat single-arm discrepancy figure adds little beyond the mathematical explanation that the marginal prior uncertainty is the same at each profile. It should justify the space it occupies.
- The application figure and table occupy float pages with substantial unused space; improve placement and sizing without compressing labels.
- A compact table defining the different borrowing diagnostics, their conditioning, units and uses would help readers follow the paper.

## Boundaries of this assessment

The review establishes specific objections through source reading, code inspection, primary-source verification, and small in-memory numerical calculations. It does not certify full implementation correctness, MCMC mixing or convergence, clinical qualification, or every numerical result. It does not recommend replacing or silently rerunning the existing results. Any change to the calibrated quantity, selection rule, or application adjustment set requires an explicit assessment of which results must be regenerated.
