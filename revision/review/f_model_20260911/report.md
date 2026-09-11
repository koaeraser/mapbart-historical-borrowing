# Explicit f model and local ESS interpretation

Section 2.2 now displays the sum-of-trees expansion for f, defines its tree, reached-leaf and leaf-mean notation, and specifies the Gaussian leaf prior and range-rule scale. It distinguishes a profile's conditional mean from an individual outcome and states the log-time interpretation for survival. The existing simulation and application tree counts are distinguished. Section 2.3 uses the same f-leaf symbol in its conditional draw.

Section 2.4 defines V_f(x) as the RWD-only posterior variance of the conditional mean. Local ESS equates the working variance, including the allowed source discrepancy, with the variance of a mean of hypothetical independent controls at the same profile. Survival reference units are uncensored log times. The single-arm compatibility limitation is retained.

No fitted model, calibration convention, numerical result, table, figure or citation changed. All existing numbered equations retain their numbers, including the leaf mixture (2). The new f display is unnumbered. Author fields remain blank.

The clean PDF has 80 pages and the cumulative tracked PDF has 101 pages. Both recorder-enabled latexmk builds complete without warnings, undefined references/citations or overfull boxes. Full-document overviews and the model, posterior-update and ESS passages were visually inspected. The cumulative redline retains the original pre-Codex baseline and the existing table-aware generator.

Integrity checker: citations, labels, sections, numbered equations, protected blocks, graphics and structure pass. Its one reference failure is the authorized extra reference to equation (1) in the conditional-mean explanation. Audited math/number additions are the explicit f expansion and prior scale, the application tree count already given in Section 5, the f-leaf symbol, and the local variance/reference definition. The full diff and checker output are retained.
