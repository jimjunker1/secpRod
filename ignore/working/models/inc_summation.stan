data {
  int<lower=1> T;              // number of time points
  int<lower=1> R;              // number of replicates per time point
  matrix[T, R] W_obs;          // replicate mean mass observations
  matrix[T, R] N_indiv;        // total individuals in each replicate
  matrix[T, R] N_obs;          // total abundance estimate per replicate
  vector[T - 1] delta_t; // days between t-1 and t
}

parameters {
  vector<lower=0>[T] W_true;   // population-level mean mass at each time
  vector[T] log_N_true;        // population-level log abundance at each time
  real<lower=0> sigma_W;       // baseline variability in W_obs
  real<lower=0> sigma_N;       // baseline variability in N_obs
}

transformed parameters{
  vector<lower = 0>[T] N_true = exp(log_N_true);
}
model {
  // Priors
  W_true ~ normal(1, 0.5);         // weakly informative prior for mean mass
  // N_true ~ lognormal(mu_N, sigma_N);        // weakly informative prior for abundance
  log_N_true ~ normal(0, 2.5);
  sigma_W ~ cauchy(0, 1);          // half-Cauchy prior
  sigma_N ~ normal(0, 0.5);

for (t in 1:T) {
  for (r in 1:R) {
    if (W_obs[t, r] != -999) {
      W_obs[t, r] ~ normal(W_true[t], sigma_W / sqrt(N_indiv[t, r]));
    }
    if (N_obs[t, r] != -999)
    N_obs[t, r] ~ lognormal(log(N_true[t]), sigma_N);

  }
}
}

generated quantities {
  real B_initial = N_true[1] * W_true[1];
  real P = B_initial;
  vector[T - 1] growth_rate;

for (t in 2:T) {
  real dW = W_true[t] - W_true[t - 1];
  real Nbar = (N_true[t] + N_true[t - 1]) / 2.0;
  P += Nbar * dW;  // or: P += Nbar * dW / delta_t[t - 1] * delta_t[t - 1]; same result
  growth_rate[t - 1] = log(W_true[t] / W_true[t -1]) / delta_t[t - 1];
}
}
