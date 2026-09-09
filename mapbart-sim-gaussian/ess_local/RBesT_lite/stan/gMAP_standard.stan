// Trimmed gMAP Stan model.
// Only the configuration used in simulate_ess.R is retained:
//   - normal/identity-link likelihood
//   - InvGamma prior on the between-trial sd tau
//   - normal random effects, single tau stratum
//   - non-centered parametrisation
//   - posterior sampling (no prior-predictive mode)
// All branching on link / tau_prior_dist / re_dist / ncp / prior_PD
// has been dropped, and so have all binomial / poisson data inputs.

data {
  // historical trials
  int<lower=1> H;

  // normal summary data
  vector[H] y;
  vector[H] y_se;

  // groups (one tau stratum)
  int<lower=1> n_groups;
  array[H] int<lower=1, upper=n_groups> group_index;

  // intercept-only design matrix (mX = 1)
  int<lower=1> mX;
  matrix[H, mX] X;

  // priors
  matrix[mX, 2] beta_prior;       // [mean, sd] per coefficient
  vector[2]     tau_prior;        // InvGamma(alpha, beta) on tau

  // location/scale guesses for the centred-NCP rescaling
  array[2] vector[mX] beta_raw_guess;
  array[2] real       tau_raw_guess;
}

parameters {
  vector[mX]       beta_raw;
  real             tau_raw;
  vector[n_groups] xi_eta;
}

transformed parameters {
  vector[mX] beta = beta_raw_guess[1] + beta_raw_guess[2] .* beta_raw;
  real       tau  = exp(tau_raw_guess[1] + tau_raw_guess[2] * tau_raw);
  vector[H]  theta;
  for (h in 1 : H)
    theta[h] = X[h] * beta + xi_eta[group_index[h]] * tau;
}

model {
  // standardised random effect (Matt trick / NCP)
  xi_eta ~ normal(0, 1);

  // regression coefficient prior
  beta ~ normal(beta_prior[ : , 1], beta_prior[ : , 2]);

  // InvGamma prior on tau and the Jacobian for the
  // log-shift-scale transform on tau_raw
  tau ~ inv_gamma(tau_prior[1], tau_prior[2]);
  target += tau_raw_guess[2] * tau_raw;

  // data likelihood
  y ~ normal(theta, y_se);
}

generated quantities {
  // intercept-only predictive draw on the (identity) response scale
  real theta_pred      = normal_rng(beta[1], tau);
  real theta_resp_pred = theta_pred;
}
