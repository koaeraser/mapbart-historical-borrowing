// Model "hybrid" (two-arm hybrid design): same as "standard" but with mu == 0.
//
//   theta_h     = tau * xi_h,        xi_h ~ N(0, 1)
//   theta_star  = tau * xi_star,     xi_star ~ N(0, 1)
//   tau         ~ InvGamma(a, b)         (same prior as default)
//   y_h         ~ N(theta_h, y_se_h^2)
//
// There is no overall mean parameter mu; the historical means
// theta_h are centred at zero with population spread tau.
// theta_resp_pred is the predictive draw from N(0, tau^2).

data {
  int<lower=1> H;
  vector[H]    y;
  vector[H]    y_se;

  int<lower=1> n_groups;
  array[H] int<lower=1, upper=n_groups> group_index;

  vector[2]    tau_prior;       // InvGamma(alpha, beta) on tau
  array[2] real tau_raw_guess;  // log-shift-scale init for tau
}

parameters {
  real             tau_raw;
  vector[n_groups] xi_eta;
}

transformed parameters {
  real      tau = exp(tau_raw_guess[1] + tau_raw_guess[2] * tau_raw);
  vector[H] theta;
  for (h in 1 : H)
    theta[h] = xi_eta[group_index[h]] * tau;   // mu = 0
}

model {
  xi_eta ~ normal(0, 1);
  tau    ~ inv_gamma(tau_prior[1], tau_prior[2]);
  target += tau_raw_guess[2] * tau_raw;        // Jacobian for tau_raw
  y      ~ normal(theta, y_se);
}

generated quantities {
  real theta_pred      = normal_rng(0, tau);
  real theta_resp_pred = theta_pred;
}
