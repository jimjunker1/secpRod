data {
  int<lower=1> T;                     // Number of unique time points
  vector[T] time;                    // Time in days

  int<lower=1> R;                    // Number of density replicate counts
  int<lower=1, upper=T> time_idx[R]; // Time index for each replicate
  int<lower=0> count[R];             // Observed replicate counts

  int<lower=1> N_mass;               // Number of mass observations
  int<lower=1, upper=T> mass_time[N_mass];  // Time index for each mass
  vector<lower=0>[N_mass] mass;      // Individual mass observations
}

parameters {
  real<lower=0> sigma_W;             // Mass observation error
  real<lower=0> alpha_W;             // GP marginal SD for mass
  real<lower=0> rho_W;               // GP length-scale for mass
  vector[T] z_W;                     // Latent GP for mass

  real<lower=0> alpha_N;             // GP marginal SD for density
  real<lower=0> rho_N;               // GP length-scale for density
  vector[T] z_N;                     // Latent GP for log-density
}

transformed parameters {
  matrix[T, T] K_W = cov_exp_quad(time, alpha_W, rho_W);
  matrix[T, T] L_W = cholesky_decompose(add_diag(K_W, 1e-6));
  vector[T] W_true = L_W * z_W;

  matrix[T, T] K_N = cov_exp_quad(time, alpha_N, rho_N);
  matrix[T, T] L_N = cholesky_decompose(add_diag(K_N, 1e-6));
  vector[T] N_true = exp(L_N * z_N);
}

model {
  // Priors
  sigma_W ~ exponential(1);
  alpha_W ~ normal(0, 2);
  rho_W ~ normal(0, 20);
  alpha_N ~ normal(0, 100);
  rho_N ~ normal(0, 20);

  z_W ~ normal(0, 1);
  z_N ~ normal(0, 1);

  // Likelihood: masses
  for (i in 1:N_mass)
    mass[i] ~ normal(W_true[mass_time[i]], sigma_W);

  // Likelihood: counts
  for (r in 1:R)
    count[r] ~ poisson(N_true[time_idx[r]]);
}

generated quantities {
  real P = 0;
  for (t in 1:(T - 1)) {
    P += 0.5 * (N_true[t] + N_true[t + 1]) * (W_true[t + 1] - W_true[t]);
  }
}
