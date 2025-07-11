library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Use the sampled summary data
obs_data <- summary_stats %>%
  select(day, mean_mass, mean_larval_density) %>%
  filter(!is.na(mean_mass), !is.na(mean_larval_density))  # optional clean-up

stan_data <- list(
  T = nrow(obs_data),
  time_obs = obs_data$day,
  W_obs = obs_data$mean_mass,
  N_obs = obs_data$mean_larval_density
)
dput(stan_data)
# Load or define the Stan code (from the parametric version)
stan_model_file <- here::here("ignore/working/models/parametric_vbg_exp_production_model.stan")
model <- stan_model(stan_model_file)
fit <- sampling(model, data = stan_data, iter = 1000, chains = 4, seed = 1312)

posterior <- extract(fit)
time <- stan_data$time_obs

# Plot fitted vs observed
W_fit <- apply(posterior$W_true, 2, mean)
N_fit <- apply(posterior$N_true, 2, mean)

plot_df <- tibble(
  day = time,
  W_obs = stan_data$W_obs,
  N_obs = stan_data$N_obs,
  W_fit = W_fit,
  N_fit = N_fit
)

dput(plot_df)
library(ggplot2)

# Mass fit
ggplot(plot_df, aes(x = day)) +
  geom_line(aes(y = W_fit), color = "red") +
  geom_point(aes(y = W_obs), color = "black") +
  labs(y = "Mean mass (mg)", title = "von Bertalanffy fit to mean mass") +
  theme_minimal()

# Density fit
ggplot(plot_df, aes(x = day)) +
  geom_line(aes(y = N_fit), color = "blue") +
  geom_point(aes(y = N_obs), color = "black") +
  labs(y = "Larval density", title = "Exponential fit to larval density") +
  theme_minimal()

print(fit, pars = c("P", "M_inf", "k", "z", "sigma_W", "sigma_N"), probs = c(0.025, 0.5, 0.975))


# Recalculate true production over all days (including unobserved time points)
true_prod_data <- all_days %>%
  filter(day >= min(obs_data$day), day <= max(obs_data$day)) %>%
  filter(alive & !adult) %>%
  group_by(day) %>%
  summarise(
    N_true = n()/(20*20),
    W_true = mean(mass),
    .groups = "drop"
  ) %>%
  arrange(day)

# Increment-summation over true population
true_production <- sum(
  0.5 * (head(true_prod_data$N_true, -1) + tail(true_prod_data$N_true, -1)) *
    (tail(true_prod_data$W_true, -1) - head(true_prod_data$W_true, -1))
)

posterior_prod <- posterior$P


tibble(P = posterior_prod) %>%
  ggplot(aes(x = P)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  geom_vline(xintercept = true_production, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Posterior distribution of secondary production",
       subtitle = "Red line = true simulated value",
       x = "Secondary production (mg per m²)",
       y = "Posterior density") +
  theme_minimal()

# Numeric comparison
quantile(posterior_prod, probs = c(0.025, 0.5, 0.975))
cat("True secondary production:", true_production, "\n")

### spline model ####


#
stan_data <- list(
  T = nrow(obs_data),
  time_obs = obs_data$day,
  W_obs = obs_data$mean_mass,
  N_obs = obs_data$mean_larval_density
)

library(splines)
# Choose number of internal knots or degrees of freedom (df)
B <- bs(stan_data$time_obs, degree = 3, df = 5, intercept = TRUE)

# Create Stan data list
stan_data <- list(
  T = nrow(obs_data),
  time_obs = obs_data$day,
  K = ncol(B),
  B = unclass(B),   # Convert to matrix
  W_obs = obs_data$mean_mass,    # your observed mean mass
  N_obs = obs_data$mean_larval_density     # your observed larval density
)

# # Create basis object using mgcv's cubic spline (P-spline, bs = "ps")
# library(mgcv)
#
# s_obj <- smoothCon(s(time_obs, bs = "ps", k = 5), data = stan_data, absorb.cons = TRUE)[[1]]
# B <- s_obj$X  # Design matrix
#
# # Stan input
# stan_data <- list(
#   T = length(time_obs),
#   K = ncol(B),
#   B = B,
#   W_obs = W_obs,
#   N_obs = N_obs
# )

stan_model_file <- here::here("ignore/working/models/spline_increment_production_model.stan")
model <- stan_model(stan_model_file)
fit <- sampling(model, data = stan_data, iter = 1000, chains = 4, seed = 1312)

posterior <- extract(fit)
time <- stan_data$time_obs

# Plot fitted vs observed
W_fit <- apply(posterior$W_true, 2, mean)
N_fit <- apply(posterior$N_true, 2, mean)

plot_df <- tibble(
  day = time,
  W_obs = stan_data$W_obs,
  N_obs = stan_data$N_obs,
  W_fit = W_fit,
  N_fit = N_fit
)

# dput(plot_df)
# library(ggplot2)

# Mass fit
plot_df %>%
ggplot(aes(x = day)) +
  geom_line(aes(y = W_fit), color = "red") +
  geom_point(aes(y = W_obs), color = "black") +
  labs(y = "Mean mass (mg)", title = "von Bertalanffy fit to mean mass") +
  theme_minimal()

# Density fit
ggplot(plot_df, aes(x = day)) +
  geom_line(aes(y = N_fit), color = "blue") +
  geom_point(aes(y = N_obs), color = "black") +
  labs(y = "Larval density", title = "Exponential fit to larval density") +
  theme_minimal()
