// =============================
// MODEL 1: Linear Growth
// =============================
data {
  int<lower=1> T;
  int<lower=1> R;
  matrix[T, R] W_obs;
  matrix[T, R] N_obs;
  matrix[T, R] N_indiv;
  vector[T - 1] delta_t;
}

parameters {
  real<lower=0> W1;
  real<lower=0> g; // daily growth rate
  vector<lower=0>[T] N_true;
  real<lower=0> sigma_W;
  real<lower=0> sigma_N;
}

transformed parameters {
  vector[T] W_true;
  W_true[1] = W1;
  for (t in 2:T)
    W_true[t] = W_true[t - 1] + g * delta_t[t - 1];
}

model {
  W1 ~ normal(0.1, 0.1);
  g ~ normal(0.001, 0.001);
  N_true ~ normal(100, 50);
  sigma_W ~ cauchy(0, 1);
  sigma_N ~ cauchy(0, 10);

  for (t in 1:T) {
    for (r in 1:R) {
      if (W_obs[t, r] != -999)
        W_obs[t, r] ~ normal(W_true[t], sigma_W / sqrt(N_indiv[t, r]));
      if (N_obs[t, r] != -999)
        N_obs[t, r] ~ normal(N_true[t], sigma_N);
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
    P += Nbar * dW;
    growth_rate[t - 1] = log(W_true[t] / W_true[t - 1]) / delta_t[t - 1];
  }
}
