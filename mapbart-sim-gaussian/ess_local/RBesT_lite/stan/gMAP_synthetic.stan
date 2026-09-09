// Model "synthetic" (single-arm synthetic control): single historical
// group treated as the anchor for the predictive draw.
//
//   theta_h     ~ N(beta_prior[1], beta_prior[2]^2)
//   theta_star  ~ N(theta_h, tau^2)
//   tau         ~ InvGamma(tau_prior[1], tau_prior[2])
//   y_h         ~ N(theta_h, y_se_h^2)
//
// h is always 1 in the simulate_ess.R use case; the model accepts
// H >= 1 generically and uses theta[1] as the anchor for theta_star.

data {
  int<lower=1> H;
  vector[H]    y;
  vector[H]    y_se;

  vector[2]    beta_prior;      // [mean, sd] for theta_h prior
  vector[2]    tau_prior;       // InvGamma(alpha, beta) on tau
  array[2] real tau_raw_guess;  // log-shift-scale init for tau
}

parameters {
  vector[H] theta;              // historical-group means
  real      tau_raw;
}

transformed parameters {
  real tau = exp(tau_raw_guess[1] + tau_raw_guess[2] * tau_raw);
}

model {
  theta ~ normal(beta_prior[1], beta_prior[2]);
  tau   ~ inv_gamma(tau_prior[1], tau_prior[2]);
  target += tau_raw_guess[2] * tau_raw;        // Jacobian for tau_raw
  y     ~ normal(theta, y_se);
}

generated quantities {
  real theta_pred      = normal_rng(theta[1], tau);
  real theta_resp_pred = theta_pred;
}
