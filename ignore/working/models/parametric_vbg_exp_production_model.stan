// Parametric model: von Bertalanffy growth and exponential mortality
// estimatees production from sampling assuming von Bertalanffy growth
// and negative exponential mortality

data {
  int<lower=2> T;                 // number of time points
  vector[T] time_obs;             // days sampled
  vector[T] W_obs;               // observed mean mass
  vector[T] N_obs;               // observed density
}

parameters {
  real<lower=0> M_inf;            // asymptotic mass
  real<lower=0> k;                // growth coefficient
  real<lower=0> t0;               // growth time offset

  real<lower=0> N0;               // initial density
  real<lower=0> z;                // mortality rate

  real<lower=0> sigma_W;          // observation error for mass
  real<lower=0> sigma_N;          // observation error for density
}

transformed parameters {
  vector[T] W_true;
  vector[T] N_true;
  for (t in 1:T) {
    W_true[t] = M_inf * (1 - exp(-k * (time_obs[t] - t0)));
    N_true[t] = N0 * exp(-z * (time_obs[t] - time_obs[1]));
  }
}

model {
  // Priors
  M_inf ~ normal(5, 2);
  k ~ normal(0.01, 0.01);
  t0 ~ normal(0, 50);
  N0 ~ normal(500, 100);
  z ~ normal(0.01, 0.01);

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
