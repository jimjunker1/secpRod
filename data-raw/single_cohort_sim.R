## code to prepare `single_cohort_sim` dataset goes here
# Spatial growth and mortality simulation in R
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggdist)

# Parameters
grid_size <- 20
mu_N_init <- 500
sigma_N_init <- 100
initial_mass <- 0.0006
mu_ln <- log(5^2 / sqrt(0.5^2 + 5^2))
sigma_ln <- sqrt(log(1 + (0.5^2 / 5^2)))
mu_z <- 0.04
sigma_z <- 0.01
cpi_start <- 290
cpi_end <- 310
days <- 506
sample_interval <- 30
sample_start <- 1    # adjustable start day
sample_end <- 365    # adjustable end day
S <- 10  # number of cells to sample per event

# Function to initialize a cohort
init_cohort <- function(i, j, start_day) {
  N_init <- max(1, round(rnorm(1, mu_N_init, sigma_N_init)))
  M_inf <- rlnorm(N_init, meanlog = mu_ln, sdlog = sigma_ln)
  k <- log(M_inf / initial_mass) / runif(N_init, cpi_start, cpi_end)
  z <- rnorm(1, mu_z, sigma_z)
  tibble(
    x = i,
    y = j,
    id = 1:N_init,
    mass = initial_mass,
    M_inf = M_inf,
    k = k,
    alive = TRUE,
    adult = FALSE,
    z = z,
    cohort_start = start_day
  )
}

# Initialize first cohort on day 1
set.seed(1312)

grid_population <- map2_dfr(rep(1:grid_size, each = grid_size), rep(1:grid_size, times = grid_size), ~init_cohort(.x, .y, 1))

grid_population %>%
  summarise(n=n(), .by = c('x','y')) %>%
  ggplot()+
  geom_tile(aes(x = x, y = y, fill = n))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  viridis::scale_fill_viridis()


# Daily update function for larval individuals only
update_day <- function(pop, current_day) {
  pop %>%
    filter(alive & !adult) %>%
    mutate(
      time_since_start = current_day - cohort_start,
      alive = runif(n()) > z,
      mass = M_inf * (1 - exp(-k * time_since_start)),
      adult = mass >= M_inf*0.98
    ) %>%
    filter(alive & !adult)  # remove those who died or became adult
}

# Run simulation with second cohort added at day 366
simulation <- vector("list", length = days)
simulation[[1]] <- grid_population

#set seed to reproduce
for (d in 2:days) {
  updated_pop <- update_day(simulation[[d - 1]], d - 1)

  if (d == 366) {
    new_cohort <- map2_dfr(rep(1:grid_size, each = grid_size), rep(1:grid_size, times = grid_size), ~init_cohort(.x, .y, 366))
    updated_pop <- bind_rows(updated_pop, new_cohort)
  }
  simulation[[d]] <- updated_pop
}

# Combine for sampling
all_days <- bind_rows(simulation, .id = "day") %>%
  mutate(day = as.integer(day))

# Sampling protocol (larvae only)
sampling_results <- list()
sampled_cells <- list()

for (t in seq(sample_start, sample_end, by = sample_interval)) {
  all_cells <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  if (length(sampled_cells) > 0) {
    prev_sampled <- bind_rows(sampled_cells)
    available_cells <- anti_join(all_cells, prev_sampled, by = c("x", "y"))
  } else {
    available_cells <- all_cells
  }
  sampled <- available_cells %>% sample_n(min(S, nrow(available_cells)))
  sampled_cells[[length(sampled_cells) + 1]] <- sampled

  sampled_data <- all_days %>%
    filter(day == t) %>%
    semi_join(sampled, by = c("x", "y")) %>%
    group_by(x, y) %>%
    summarise(
      larvalDensity = n(),
      massDistribution = list(mass),
      .groups = "drop"
    ) %>%
    mutate(day = t)

  sampling_results[[length(sampling_results) + 1]] <- sampled_data
}

# Final output
daily_sampling <- bind_rows(sampling_results)
singleCohortSim <- daily_sampling

usethis::use_data(singleCohortSim, overwrite = TRUE)


### Below was used to summarise and visualize the data sets ####

# load(here::here("ignore/working/sampling and model test/single_cohort_sim.RData"))
# # Summarize samples over time
# summary_stats <- daily_sampling %>%
#   unnest(massDistribution) %>%
#   group_by(day) %>%
#   summarise(
#     mean_mass = mean(massDistribution, na.rm = TRUE),
#     sd_mass = sd(massDistribution, na.rm = TRUE),
#     mean_larvalDensity = mean(larvalDensity),
#     .groups = "drop"
#   )
#
# # Plot mean larval density and mass
# ggplot(summary_stats, aes(x = day)) +
#   stat_halfeye(data = daily_sampling, aes(x = day, y = larvalDensity),
#                color = 'green')+
#   stat_halfeye(data = daily_sampling %>% unnest(massDistribution), aes(x = day, y = massDistribution*100),
#                color = 'red')+
#   geom_line(aes(y = mean_larvalDensity), color = 'green') +
#   geom_line(aes(y = mean_mass * 100), color = 'red') +
#   scale_y_continuous(
#     name = "Larval Density",
#     sec.axis = sec_axis(~./100, name = "Mean Mass (mg)")
#   ) +
#   theme_minimal() +
#   labs(title = "Larval Density and Mean Mass over Time", x = "Day")
#
# # Ridge plot of larval mass distributions over time
# daily_sampling %>%
#   unnest(massDistribution) %>%
#   ggplot(aes(x = massDistribution, y = factor(day))) +
#   stat_halfeye(.width = 0.5, fill = 'gray70') +
#   theme_minimal() +
#   labs(x = "Larval Mass (mg)", y = "Day", title = "Larval Mass Distributions Over Time")+
#   coord_flip()
#

