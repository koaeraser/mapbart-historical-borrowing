# Phase 3 critic report: 03-theory/appendix_new.tex

Date: 2026-09-09. Every Lemma, Proposition, Theorem and Corollary checked step by step against method_spec.md (A1, B, C, D) and, where relevant, 01-code/lrcbart. Constants recomputed in R.

## Summary table

| Result | Claim | Verdict | Reason |
|---|---|---|---|
| Lemma A1 | Exact f and g leaf marginals and conditionals | CORRECT | M_f, M_g, Bernoulli weight and conditionals verified by completing the square; matches spec B(i), B(iii). |
| Proposition A1 | Tree moves plus leaf draws target the exact posterior | MINOR | Argument right; sampler truncates tau_0^2 to (0, tau_1^2), so the prior must be stated as truncated. |
| Theorem A1 (T1) | Bounded RWD-induced shift, Gaussian decay | CORRECT | logit p = a + kappa ybar^2 verified; sup lemma verified in both cases; 0.51, 0.31 (at 0.43), 0.54 recomputed as 0.509, 0.306 (at 0.434), 0.540. |
| Corollary A1 | Bounded shift of the standardized control mean | CORRECT | Woodbury step, mean-difference norm tau_1/tau_0, determinant ratio, union over 2^{L_g}-1 configurations verified; (c) bound 1.80 and global exact sup 0.067 at 0.30 reproduced. |
| Theorem A2 (T2) | Leaf detection, half-detection point, large-n limits | CORRECT | d_{1/2} = 0.540; abstention 0.011, 0.43, 0.98, 1.00; limit at d = 0 is 0.993. Gap paragraph correctly retracts the spec's false limit claim. |
| Theorem A3 | Discrepancy map, single tree | CORRECT | Mixture form, both bounds and consistency verified. |
| Proposition A2 | Discrepancy map, H_g trees, fixed partitions | CORRECT | Sigma_z dominated by D_z, cell-mean variance, Cauchy-Schwarz bias and noise bounds verified. |
| Theorem A4(a) (T3) | ESS_0 continuous, convex, strictly decreasing, limits | MINOR | Untruncated case correct (interchange justified by 1/(4Vs)); "if truncated, the same conclusions hold" is false for the s to infinity limit. |
| Theorem A4(b) | Half-line sublevel sets, G a CDF, feasibility | CORRECT | P(N_target < sigma_1^2/V) > 1 - q is exactly s_star > 0. |
| Theorem A4(c) | Monte Carlo consistency | MINOR | Needs aperiodicity; block length m unstated; implementation runs fixed B = 20 with m growing, the opposite regime. |
| Theorem A4(d) | Anchor n_2 sigma_1^2/sigma_2^2 | CORRECT | Cauchy-Schwarz and equality case verified. |
| Proposition A3 | Population identification of f and g | CORRECT | Consistent with the model on all four support regions; single-arm variance correct. |

Zero MAJOR, three MINOR.

## MINOR items and fixes

**Proposition A1, prior of tau_0^2.** Line 46 says the sampler draws tau_0^2 truncated to (0, tau_1^2); `update_tau0_cpp` confirms it (`rsinvchisq_upper`). Equation (leafprior) and spec A1 give an untruncated scaled-Inv-chi^2. A truncated conditional is conjugate for the truncated prior only, so the "exact joint posterior" is the one under the truncated prior. Fix: state the prior as truncated in Web Appendix A; keep the Appendix D sentence that calibration uses the untruncated prior, adding that the truncation mass is negligible at s_0^2 of order 10^{-3} against tau_1^2 of order 1.

**Theorem A4(a), truncated prior.** Under truncation, as s to infinity the density, proportional to t^{-nu_0/2-1} exp(-nu_0 s/(2t)) on (0, tau_1^2), concentrates at tau_1^2, so ESS_0 tends to sigma_1^2/(V + c tau_1^2) > 0, not 0; the range is (sigma_1^2/(V + c tau_1^2), sigma_1^2/V). Convexity and differentiability rest on tau_0^2 = nu_0 s/Q, which fails under truncation. Fix: replace "the same conclusions hold" by "ESS_0 remains continuous and strictly decreasing, with limit sigma_1^2/(V + c_g H_g tau_1^2) as s to infinity"; note in (b) that under truncation a root exists only if N_target also exceeds that floor. The monotone-likelihood-ratio argument as written suffices for this.

**Theorem A4(c), block statistics.** (i) A stationary ergodic chain can be periodic, and then the m-block sequence is not ergodic (deterministic two-state alternation, m = 2); add "aperiodic". (ii) The law P of the block statistic, hence G and s_star, depends on m; write P_m. (iii) `lrc_ess_calibrate` fixes `n_blocks = 20` with block length nd/20, so the implementation's asymptotics are m to infinity at fixed B. There each V^{(b)} tends to V_mu^f, G becomes a step at s^dagger(V_mu^f), and the Pr-rule reduces to the root of ESS_0(V_mu^f; s) = N_target; the rescaling by V_whole/mean(V_b_raw) is a finite-m correction matching the block mean to the whole-chain variance and nothing more. The rescaling argument given is correct for the B to infinity regime; the CLT is not "verbatim" because gamma_B is random and needs one delta-method sentence. Fix: state the regime used, add the one-line m to infinity limit (within-block variances of a geometrically ergodic chain converge to the stationary variance), keep B to infinity as the alternative.

## Editorial (no verdict change)

Theorem A1(b)-(d) and Corollary A1(a) need 0 < w < 1 stated. Theorem A2(b): at a = 0, P_spike equals 1/2 at ybar = 0, so "less than" should be "at most". Corollary A1(a): ||a||_2^2 <= max_j a_j sum_j a_j <= H_g gives ||a||_2 <= sqrt(H_g) for free. Remarks A4 and A6 use H_g = 10 while spec A3 fixes H_g = 5; at H_g = 5 the common-partition bound is 1.15 per leaf and the exact supremum 0.148 at ybar = 0.33. The index convention (j = 1 slab in rho_j, z = 1 spike) deserves one sentence. The Gap paragraphs are honest and no statement claims more than is proved.

## Verdict

THEORY GATE PASSED
