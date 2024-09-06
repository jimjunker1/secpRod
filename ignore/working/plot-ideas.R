# Load required libraries
library(ggplot2)
library(dplyr)

# Sample data frame (replace this with your actual data)
set.seed(42)
data <- data.frame(
  date = rep(seq.Date(as.Date("2023-01-01"), by = "month", length.out = 12), each = 50),
  species = rep(c("Hydropsyche_elissoma", "Hydropsyche_incommoda", "Macronema_carolina"), times = 200),
  count = c(rnorm(200, mean = 50, sd = 10), rnorm(200, mean = 30, sd = 8), rnorm(200, mean = 20, sd = 5))
)

# Bin the data for histograms
data_binned <- data %>%
  group_by(date, species) %>%
  mutate(bin = cut(count, breaks = 30)) %>%
  count(date, species, bin) %>%
  mutate(midpoint = as.numeric(sub("\\((.+),.*", "\\1", bin)) + 5)  # Estimate midpoint of each bin

# Create the step-like plot with density curves
ggplot(data_binned, aes(x = midpoint, y = n, color = species)) +
  geom_step(direction = "hv") +  # Use horizontal-vertical steps
  facet_grid(species ~ date, scales = "free_y", space = "free_y") +  # Separate by species and date
  theme_minimal() +
  labs(x = "Count", y = "Density", title = "Stepwise Density Progression of Species over Time") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank()
  )
