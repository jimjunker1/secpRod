data {
  int<lower=1> T;                     // Number of unique time points
  int<lower=1> R;                     // Number of replicate observations
  int<lower=1, upper=T> time_idx[R];  // Time index for each replicate count
  int<lower=0> count[R];              // Observed count (density proxy)

  int<lower=1> N_mass;                // Number of observed masses
  int<lower=1, upper=T> mass_time[N_mass]; // Time index for each mass
  vector<lower=0>[N_mass] mass;       // Observed individual masses

  vector[T] time;                     // Time in days (e.g., 1:365)
}

parameters {
  real<lower=0> M_inf;                // Asymptotic mass
  real<lower=0> k;                    // Growth coefficient
  real<lower=0> N0;                   // Initial density
  real<lower=0> z;                    // Mortality rate

  real<lower=0> sigma_W;              // Mass observation SD
  real<lower=0> sigma_N;              // Density observation SD
}

transformed parameters {
  vector[T] W_true;
  vector[T] N_true;

  for (t in 1:T) {
    W_true[t] = M_inf * (1 - exp(-k * time[t]));
    N_true[t] = N0 * exp(-z * time[t]);
  }
}

model {
  // Priors
  M_inf ~ lognormal(log(5), 0.3);
  k ~ lognormal(log(0.02), 0.5);
  N0 ~ normal(500, 100);
  z ~ lognormal(log(0.01), 0.3);

  sigma_W ~ exponential(1);
  sigma_N ~ exponential(1);

  // Mass likelihood
  for (i in 1:N_mass) {
    mass[i] ~ normal(W_true[mass_time[i]], sigma_W);
  }

  // Density likelihood
  for (r in 1:R) {
    count[r] ~ poisson_log(log(N_true[time_idx[r]]));
  }
}

generated quantities {
  real P = 0;
  for (t in 1:(T - 1)) {
    real dW = W_true[t + 1] - W_true[t];
    real mean_N = 0.5 * (N_true[t] + N_true[t + 1]);
    P += mean_N * dW;
  }
}
