data {
  int<lower=2> T;                  // number of time points
  vector[T] time;                 // observed time points (e.g., day)
  vector[T] W_obs;                // observed mean mass
  vector[T] N_obs;                // observed larval density
}

parameters {
  real<lower=0> sigma_W;          // obs error: mass
  real<lower=0> sigma_N;          // obs error: density

  real<lower=0> alpha_W;          // GP marginal SD: mass
  real<lower=0> rho_W;            // GP length scale: mass

  real<lower=0> alpha_N;          // GP marginal SD: density
  real<lower=0> rho_N;            // GP length scale: density

  vector[T] z_W;                  // latent GP for mass
  vector[T] z_N;                  // latent GP for density
}

transformed parameters {
  matrix[T, T] K_W = cov_exp_quad(time, alpha_W, rho_W);
  matrix[T, T] L_W = cholesky_decompose(add_diag(K_W, 1e-6));

  matrix[T, T] K_N = cov_exp_quad(time, alpha_N, rho_N);
  matrix[T, T] L_N = cholesky_decompose(add_diag(K_N, 1e-6));

  vector[T] W_true = L_W * z_W;
  vector[T] N_true = L_N * z_N;
}

model {
  // Priors
  alpha_W ~ normal(0, 1);
  rho_W ~ normal(0, 20);

  alpha_N ~ normal(0, 100);
  rho_N ~ normal(0, 20);

  sigma_W ~ exponential(1);
  sigma_N ~ exponential(1);

  z_W ~ normal(0, 1);
  z_N ~ normal(0, 1);

  // Observation model
  W_obs ~ normal(W_true, sigma_W);
  N_obs ~ normal(N_true, sigma_N);
}

generated quantities {
  real P = 0;
  for (t in 1:(T - 1)) {
    P += 0.5 * (N_true[t] + N_true[t + 1]) * (W_true[t + 1] - W_true[t]);
  }
}
