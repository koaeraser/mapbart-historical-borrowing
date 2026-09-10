# Primary-source and code verification

- Neuenschwander et al., Predictively Consistent Prior Effective Sample Sizes, Biometrics 76 (2020), 578-587; primary author preprint https://arxiv.org/pdf/1907.04185, Section 2.4 equation (7) and Section 2.5.1. ELIR averages the curvature of the log marginal prior divided by unit likelihood information. Normal prior yields variance ratio; a scale mixture generally does not. The direct mixture check in checks/mathematical_checks.R independently demonstrates the distinction.
- Jin et al., Bayesian adaptive design for covariate-adaptive historical control information borrowing, Statistics in Medicine (2023); primary author PDF https://jinhuaqing.github.io/files/2023_SIM_CAHB.pdf, pages 5341 and 5343. Their outcome/precision function estimates are nonparametric kernel estimates. This corrects the paper's characterization as parametric.
- lrcbart.R:323-324 and 434-437 implement an untruncated scale draw with a plug-in mean tree-sharing coefficient. The whole-chain cap at 426-443 bypasses the block Pr-rule.
- lrcbart.R:469-478 contains the tree-binomial prior helper; no calls occur in the result-generating R scripts. It is not exact for independent leaf states of a standardized estimand.
- lrcbart.R:500-508 and lrcbart.cpp:759-788 compute the reported realized diagnostic using sum_h sum_l a_hl^2 tau_z^2; the diagnostic uses the prior covariance evaluated at posterior states, not the posterior variance of g.
- 01-code/R/fit_lrcbart.R:110-120 exports that exact leaf-weighted diagnostic to the study summaries. 02-validation/sc4_boundary/pointwise_ess.R:105-113 uses the same method for regional ESS.
- 04-application/code/run_lrcbart.R:65-69 uses that diagnostic and the pointwise function with mean(tau0_sq). lrcbart.R:514-518 places the mean scale outside inversion for the application map; the extra Sc4 pointwise prior-ESS analysis instead averages over untruncated prior scale draws.
- 04-application/code/prep_data.R:15 and 40-41 uses transplant_off_protocol as the ASCT covariate. Its timing is after induction as reported in the current manuscript; the resulting analyses are descriptive adjusted associations, not baseline-adjusted causal total effects. No data/model were changed.

Source verification does not certify the full model, clinical study, sampler mixing or confirmatory operating characteristics.
