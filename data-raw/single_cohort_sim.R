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
mu_z <- 0.035
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

# grid_population %>%
#   summarise(n=n(), .by = c('x','y')) %>%
#   ggplot()+
#   geom_tile(aes(x = x, y = y, fill = n))+
#   scale_x_continuous(expand = c(0,0))+
#   scale_y_continuous(expand = c(0,0))+
#   viridis::scale_fill_viridis()


# Daily update function for larval individuals only
update_day <- function(pop, current_day) {
  pop %>%
    dplyr::filter(alive & !adult) %>%
    dplyr::mutate(
      time_since_start = current_day - cohort_start,
      alive = runif(n()) > z,
      mass = M_inf * (1 - exp(-k * time_since_start)),
      adult = mass >= M_inf*0.9999
    ) %>%
    dplyr::filter(alive & !adult)  # remove those who died or became adult
}

# Run simulation with second cohort added at day 366
simulation <- vector("list", length = days)
simulation[[1]] <- grid_population

#set seed to reproduce
for (d in 2:days) {
  updated_pop <- update_day(simulation[[d - 1]], d - 1)

  # if (d == 366) {
  #   new_cohort <- map2_dfr(rep(1:grid_size, each = grid_size), rep(1:grid_size, times = grid_size), ~init_cohort(.x, .y, 366))
  #   updated_pop <- bind_rows(updated_pop, new_cohort)
  # }
  simulation[[d]] <- updated_pop
}

# Combine for sampling
all_days <- bind_rows(simulation, .id = "day") %>%
  dplyr::mutate(day = as.integer(day))

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
    dplyr::filter(day == t) %>%
    semi_join(sampled, by = c("x", "y")) %>%
    group_by(x, y) %>%
    summarise(
      larvalDensity = n(),
      massDistribution = list(mass),
      .groups = "drop"
    ) %>%
    dplyr::mutate(day = t)

  sampling_results[[length(sampling_results) + 1]] <- sampled_data
}

# Final output
daily_sampling <- bind_rows(sampling_results) %>%
  dplyr::mutate(taxonID = 'sppX', repID = 1:n(), .by = c('day')) %>%
  rename(density = larvalDensity, afdm_mg = massDistribution, dateID = day)

zeroFills = expand.grid(taxonID = 'sppX',
                        dateID = unique(daily_sampling$dateID),
                        repID = 1:10,
                        density = 0,
                        afdm_mg = NA)

daily_sampling = daily_sampling %>%
  select(-x,-y) %>%
  full_join(.,zeroFills, by = c('taxonID','dateID','repID')) %>%
  mutate(density = ifelse(is.na(density.x), density.y, density.x),
         afdm_mg = ifelse(is.na(afdm_mg.x), afdm_mg.y, afdm_mg.x)) %>%
  select(taxonID, dateID, repID, density, afdm_mg)

singleCohortSim <- daily_sampling

usethis::use_data(singleCohortSim, overwrite = TRUE)


### Below was used to summarise and visualize the data sets ####

# # Summarize samples over time
# summary_stats <- daily_sampling %>%
#   unnest(afdm_mg, keep_empty= TRUE) %>%
#   group_by(dateID) %>%
#   summarise(
#     meanMass = mean(afdm_mg, na.rm = TRUE),
#     sdMass = sd(afdm_mg, na.rm = TRUE),
#     meanDensity = mean(density, na.rm = TRUE)
#   )
#
# # Plot mean larval density and mass
# ggplot(summary_stats, aes(x = dateID)) +
#   stat_halfeye(data = daily_sampling, aes(x = dateID, y = density),
#                color = 'green')+
#   stat_halfeye(data = daily_sampling %>% unnest(afdm_mg), aes(x = dateID, y = afdm_mg*100),
#                color = 'red')+
#   geom_line(aes(y = meanDensity), color = 'green') +
#   geom_line(aes(y = meanMass * 100), color = 'red') +
#   scale_y_continuous(
#     name = "Larval Density",
#     sec.axis = sec_axis(~./100, name = "Mean Mass (mg)")
#   ) +
#   theme_minimal() +
#   labs(title = "Larval Density and Mean Mass over Time", x = "Day")
#
# # Ridge plot of larval mass distributions over time
# daily_sampling %>%
#   unnest(afdm_mg) %>%
#   ggplot(aes(x = afdm_mg, y = factor(dateID))) +
#   stat_halfeye(.width = 0.5, fill = 'gray70') +
#   theme_minimal() +
#   labs(x = "Larval Mass (mg)", y = "Day", title = "Larval Mass Distributions Over Time")+
#   coord_flip()
