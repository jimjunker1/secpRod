# Load necessary libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(42)

# Parameters
n_individuals <- 5000
initial_mass <- 0.0006
W_inf <- 2.0
k_mean <- 0.02
k_sd <- 0.005
t0_mean <- 300 - 6/k_mean
t0_sd <- 2
mortality_rate_mean <- 2.5e-5
mortality_sd <- 0.15
days <- 0:365
dates <- as.Date("2025-01-01") + days

# Individual-specific parameters
individual_k <- rnorm(n_individuals, k_mean, k_sd)
individual_t0 <- rnorm(n_individuals, t0_mean, t0_sd)
mortality_rate <- rlnorm(n_individuals, log(mortality_rate_mean), mortality_sd)

# Simulate growth and survival
# Pre-allocate matrices
mass_matrix <- matrix(NA, nrow = n_individuals, ncol = length(days))
status_matrix <- matrix("juvenile", nrow = n_individuals, ncol = length(days))

for (i in 1:n_individuals) {
  t <- days

  # --- GROWTH ---
  growth <- W_inf * (1 - exp(-individual_k[i] * (t - individual_t0[i])))^3
  growth[t < individual_t0[i]] <- NA
  growth[growth < 0.0001 & !is.na(growth)] <- 0.0001
  growth[growth > W_inf] <- W_inf

  # --- TYPE III SURVIVAL: sample death time from CDF ---
  U <- runif(1)
  death_time <- sqrt(-log(U) / mortality_rate[i])
  death_day <- which(t >= death_time)[1]

  # Initialize as juvenile
  status_matrix[i, ] <- "juvenile"

  # --- APPLY DEATH ---
  if (!is.na(death_day)) {
    status_matrix[i, death_day:length(t)] <- "dead"
  }

  # --- APPLY ADULT TRANSITION ---
  is_adult <- growth >= 0.99 * W_inf
  for (d in seq_along(t)) {
    if (status_matrix[i, d] == "juvenile" &&
        !is.na(is_adult[d]) && is_adult[d]) {
      status_matrix[i, d:length(t)] <- "adult"
      break
    }
  }

  # --- ASSIGN MASS IF JUVENILE ONLY ---
  for (d in seq_along(t)) {
    if (status_matrix[i, d] == "juvenile") {
      mass_matrix[i, d] <- growth[d]
    }
  }
}

### Sampling simulation ----
sample_data <- list()
sample_area <- 1   # area per replicate
replicates_per_day <- 5
sample_days <- round(seq(0, 365, length.out = 11))

for (t in sample_days) {
  juvenile_ids <- which(status_matrix[, t + 1] == "juvenile")
  n_juveniles <- length(juvenile_ids)

  for (rep in 1:replicates_per_day) {
    if (n_juveniles < 10) next  # skip replicate if too few juveniles

    # Sample proportionally to system-wide juvenile density
    expected_sample_size <- rpois(1, lambda = n_juveniles / (n_individuals / sample_area))
    sample_ids <- sample(juvenile_ids, size = min(expected_sample_size, n_juveniles), replace = FALSE)

    sampled_masses <- mass_matrix[sample_ids, t + 1]
    sampled_masses <- sampled_masses[!is.na(sampled_masses)]
    if (length(sampled_masses) == 0) next

    df <- data.frame(
      taxonID = "species_1",
      repID = paste0("rep_", rep),
      dateID = format(dates[t + 1], "%Y-%m-%d"),
      massClass = round(sampled_masses, 4),
      count = 1  # each individual counts as one
    )
    sample_data[[length(sample_data) + 1]] <- df
  }
}

adult_count_by_day300 <- sum(status_matrix[, 301] == "adult")
print(adult_count_by_day300)

# Combine and aggregate to final long format
df <- bind_rows(sample_data) %>%
  group_by(taxonID, repID, dateID, massClass) %>%
  summarise(count = sum(count), .groups = "drop")

# Compute mean mass of surviving individuals over time
# Plot mean mass of surviving individuals
mean_mass_over_time <- colMeans(mass_matrix, na.rm = TRUE)

mass_df <- data.frame(
  date = as.Date("2025-01-01") + 0:365,
  mean_mass = mean_mass_over_time
)

ggplot(mass_df, aes(x = date, y = mean_mass)) +
  geom_line(color = "blue") +
  labs(title = "Mean Mass of Cohort Over Time", y = "Mean Mass (mg)", x = "Date") +
  theme_minimal()

# Plot percent alive over time
percent_alive <- colSums(!is.na(mass_matrix)) / n_individuals * 100

alive_df <- data.frame(
  date = as.Date("2025-01-01") + 0:365,
  percent_alive = percent_alive
)

ggplot(alive_df, aes(x = date, y = percent_alive)) +
  geom_line(color = "red") +
  labs(title = "Cohort Survival Over Time", y = "Percent Alive (%)", x = "Date") +
  theme_minimal()

density_df <- df %>%
  group_by(dateID) %>%
  summarise(density = sum(count), .groups = "drop") %>%
  mutate(date = as.Date(dateID))

ggplot(density_df, aes(x = date, y = density)) +
  geom_line(size = 1.2, color = "darkgreen") +
  labs(title = "Estimated Cohort Density Over Time",
       x = "Date",
       y = "Total Sampled Density (counts/unit area)") +
  theme_minimal()

true_density <- colSums(!is.na(mass_matrix))  # number of alive individuals per day
plot(as.Date("2025-01-01") + 0:365, true_density, type = "l")

# True juvenile density
juvenile_density <- colSums(status_matrix == "juvenile")
adult_cumulative <- colSums(status_matrix == "adult")

date_seq <- as.Date("2025-01-01") + 0:365

# Plot juvenile density
plot_df1 <- data.frame(date = date_seq, juveniles = juvenile_density)
ggplot(plot_df1, aes(x = date, y = juveniles)) +
  geom_line(color = "darkgreen") +
  labs(title = "Juvenile Density Over Time", y = "Individuals", x = "Date") +
  theme_minimal()

# Plot cumulative adult emergence
plot_df2 <- data.frame(date = date_seq, adults = adult_cumulative)
ggplot(plot_df2, aes(x = date, y = adults)) +
  geom_line(color = "purple") +
  labs(title = "Cumulative Adult Emergence", y = "Individuals", x = "Date") +
  theme_minimal()

# Summarise for Stan input
df_summary <- df %>%
  group_by(dateID, repID) %>%
  summarise(
    W_obs = sum(massClass * count) / sum(count),
    N_obs = sum(count),
    N_indiv = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    date = as.integer(as.factor(dateID)),
    rep = as.integer(as.factor(repID))
  )

T <- length(unique(df_summary$date))
R <- max(table(df_summary$date))

# Create Stan matrices
W_mat <- matrix(NA, nrow = T, ncol = R)
N_mat <- matrix(NA, nrow = T, ncol = R)
N_indiv_mat <- matrix(NA, nrow = T, ncol = R)

for (t in 1:T) {
  row_t <- df_summary %>% filter(date == t)
  W_mat[t, 1:nrow(row_t)] <- row_t$W_obs
  N_mat[t, 1:nrow(row_t)] <- row_t$N_obs
  N_indiv_mat[t, 1:nrow(row_t)] <- row_t$N_indiv
}

# Compute delta_t (non-uniform spacing)
date_lookup <- df_summary %>%
  distinct(dateID, date) %>%
  arrange(date)

delta_t <- diff(as.numeric(as.Date(date_lookup$dateID)))

# Final list for Stan
stan_data <- list(
  T = T,
  R = R,
  W_obs = W_mat,
  N_obs = N_mat,
  N_indiv = N_indiv_mat,
  delta_t = delta_t
)
sentinel <- -999
stan_data$W_obs[is.na(stan_data$W_obs)] <- sentinel
stan_data$N_obs[is.na(stan_data$N_obs)] <- sentinel
stan_data$N_indiv[is.na(stan_data$N_indiv)] <- sentinel

fit <- stan(file = here("ignore/working/models/inc_summation.stan"), data = stan_data,
            chains = 4, iter = 2000, seed = 1312)

print(fit, pars = c("B_initial", "P"))
