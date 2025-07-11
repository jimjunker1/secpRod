
cohort_to_pseudotime <- function(time_vec, cohort_bounds) {
  # cohort_bounds is a list of numeric vectors giving time ranges, e.g.
  # list(c(90, 150), c(365, 500))

  pseudo_time <- numeric(length(time_vec))
  cohort_labels <- numeric(length(time_vec))

  for (i in seq_along(cohort_bounds)) {
    bounds <- cohort_bounds[[i]]
    in_cohort <- time_vec >= bounds[1] & time_vec <= bounds[2]

    # Map to pseudo-time within cohort
    pseudo_time[in_cohort] <- time_vec[in_cohort] - bounds[1] + 1
    cohort_labels[in_cohort] <- i
  }

  tibble(
    original_time = time_vec,
    pseudo_time = pseudo_time,
    cohort = cohort_labels
  )
}


revert_pseudotime <- function(pseudo_df, cohort_bounds) {
  # pseudo_df must include 'pseudo_time' and 'cohort' columns
  stopifnot(all(c("pseudo_time", "cohort") %in% names(pseudo_df)))

  pseudo_df %>%
    rowwise() %>%
    mutate(
      original_time = cohort_bounds[[cohort]][1] + pseudo_time - 1
    ) %>%
    ungroup()
}


time_obs <- c(100, 130, 370, 400, 430, 460)
cohort_bounds <- list(c(90, 150), c(365, 500))

pseudo_map <- cohort_to_pseudotime(time_obs, cohort_bounds)

# Add to your modeling dataframe
model_df <- bind_cols(pseudo_map, tibble(W_obs, N_obs))

# Fit model using model_df$pseudo_time and generate spline basis on that
B <- bs(model_df$pseudo_time, degree = 3, df = 5, intercept = TRUE)


model_output <- tibble(
  pseudo_time = 1:60,
  cohort = c(rep(1, 30), rep(2, 30)),
  W_fit = ..., N_fit = ...
)

# Back to original time
model_output_original <- revert_pseudotime(model_output, cohort_bounds)
