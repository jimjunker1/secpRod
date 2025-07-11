// Spline-based model: flexible smooth trajectories for mass and density using basis expansion

data {
  int<lower=2> T;                      // number of observed time points
  int<lower=1> K;                      // number of spline basis functions
  matrix[T, K] B;                      // spline basis matrix for time_obs
  vector[T] W_obs;                     // observed mean mass
  vector[T] N_obs;                     // observed larval density
}

parameters {
  vector[K] beta_W;                    // spline coefficients for mass
  vector[K] beta_N;                    // spline coefficients for density

  real<lower=0> sigma_W;               // observation error for mass
  real<lower=0> sigma_N;               // observation error for density
}

transformed parameters {
  vector[T] W_true = B * beta_W;
  vector[T] N_true = B * beta_N;
}

model {
  // Priors
  beta_W ~ normal(0, 2);
  beta_N ~ normal(0, 100);

  sigma_W ~ normal(0, 1);
  sigma_N ~ normal(0, 1);

  // Likelihood
  W_obs ~ normal(W_true, sigma_W);
  N_obs ~ normal(N_true, sigma_N);
}

generated quantities {
  real P = 0;
  for (t in 1:(T - 1)) {
    P += 0.5 * (N_true[t] + N_true[t + 1]) * (W_true[t + 1] - W_true[t]);
  }
}
