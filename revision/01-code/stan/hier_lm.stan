// hier_lm.stan: HMC version of HierLM (paper Section "Hierarchical
// Regression Models") used only to cross-check the conjugate Gibbs sampler
// in R/comparators.R (gibbs_hier). Gaussian outcome, two groups.
data {
  int<lower=1> N;
  int<lower=1> K;                 // p + 1 (intercept first)
  matrix[N, K] X;                 // design with a leading column of ones
  vector[N] y;
  array[N] int<lower=1, upper=2> g;
  vector<lower=0>[2] lambda;      // BART-style residual scale per group
  real<lower=0> nu;
}
parameters {
  array[2] vector[K] theta;
  vector[K] mu;
  vector<lower=0>[K] tau2;
  vector<lower=0>[2] sigma2;
}
model {
  mu ~ normal(0, 10);
  tau2 ~ inv_gamma(1.5, 0.75);
  sigma2 ~ inv_gamma(nu / 2, nu * lambda / 2);
  for (j in 1:2) theta[j] ~ normal(mu, sqrt(tau2));
  for (i in 1:N) y[i] ~ normal(X[i] * theta[g[i]], sqrt(sigma2[g[i]]));
}
