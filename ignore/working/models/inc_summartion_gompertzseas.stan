// =============================
// MODEL 3: Gompertz Growth with time-varying B(t) using splines
// =============================
data {
  int<lower=1> T;
  int<lower=1> R;
  int<lower=1> J;
  matrix[T, R] W_obs;
  matrix[T, R] N_obs;
  matrix[T, R] N_indiv;
  vector[T] time;
  matrix[T, J] B_spline;
}

parameters {
  real<lower=0> A;
  real<lower=0> B0;
  vector[J] beta_B;
  real<lower=0> t0;
  vector<lower=0>[T] N_true;
  real<lower=0> sigma_W;
  real<lower=0> sigma_N;
}

transformed parameters {
  vector[T] B_t;
  vector[T] W_true;
  B_t = B0 + B_spline * beta_B;
  for (t in 1:T)
    W_true[t] = A * exp(-exp(-B_t[t] * (time[t] - t0)));
}

model {
  A ~ normal(3, 1);
  B0 ~ normal(0.05, 0.02);
  beta_B ~ normal(0, 0.1);
  t0 ~ normal(150, 20);
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
  }

  for (t in 1:(T - 1)) {
    growth_rate[t] = B_t[t] * exp(-B_t[t] * (time[t] - t0));
  }
}
