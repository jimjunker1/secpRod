data {
  int<lower=1> T;                      // number of unique time points
  int<lower=1> K;                      // number of spline basis functions
  matrix[T, K] B;                      // spline basis matrix (e.g. from splines::bs)

  int<lower=1> R;                      // number of density replicate counts
  int<lower=1, upper=T> time_idx[R];   // time index for each replicate
  int<lower=0> count[R];               // observed replicate count

  int<lower=1> N_mass;                 // number of individual mass observations
  int<lower=1, upper=T> mass_time[N_mass];  // time index for each mass
  vector<lower=0>[N_mass] mass;        // individual mass observations
}

parameters {
  vector[K] beta_W;                   // spline coefficients for mass
  vector[K] beta_N;                   // spline coefficients for density

  real<lower=0> sigma_W;              // observation SD for mass
}

transformed parameters {
  vector[T] W_true = B * beta_W;      // predicted mean mass at time t
  vector[T] N_true = exp(B * beta_N); // predicted mean density at time t
}

model {
  // Priors
  beta_W ~ normal(0, 2);
  beta_N ~ normal(0, 10);

  sigma_W ~ exponential(1);

  // Likelihood: individual masses
  for (i in 1:N_mass)
    mass[i] ~ normal(W_true[mass_time[i]], sigma_W);

  // Likelihood: density replicates
  for (r in 1:R)
    count[r] ~ poisson(N_true[time_idx[r]]);
}

generated quantities {
  real P = 0;
  for (t in 1:(T - 1)) {
    P += 0.5 * (N_true[t] + N_true[t + 1]) * (W_true[t + 1] - W_true[t]);
  }
}
