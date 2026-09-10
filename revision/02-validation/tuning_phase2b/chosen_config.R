# Candidates carried into the confirmation stage after the OFAT, combo and
# moderate stages (summary_ofat.csv, summary_combo.csv, summary_moderate.csv),
# and the final choice.
# H3 with the (0.95, 2) split prior was the first candidate (best Sc4 map) and
# was dropped in the confirmation run: Sc3 bias 0.10 and map 0.17 (a global
# shift partly hidden in the spike by deeper trees), Sc4 delta 1 bias -0.085.
# Fewer trees (H1, H3) trade the Sc4 delta 2 map for Sc3 / Sc4 delta 1 bias.
CANDIDATES <- list(
  list(H_g = 3, alpha_g = 0.95, beta_g = 2),
  list(H_g = 1), list(H_g = 3), list(H_g = 5),
  list(H_g = 3, tau1_mult = 4),
  list(H_g = 10, g_sweeps = 5),
  list(H_g = 5, w_prior = c(1, 1)),
  list(H_g = 10, w_prior = c(1, 1)),
  list(H_g = 5, tau1_mult = 0.5, w_prior = c(1, 1))
)
# Final choice (set after the confirmation stage; see NOTES.md "Phase 2b tuning").
CHOSEN <- list(H_g = 5, alpha_g = 0.5, beta_g = 3, tau1_mult = 1, n_min_g = 5, g_sweeps = 1, w_prior = c(1, 1))
