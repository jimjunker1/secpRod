# Sample data for community secondary production analysis

A collection of community macroinvertebrate data for estimating
secondary production

## Usage

``` r
overlapCohortSim
```

## Format

\`overlapCohortSim\` A single tibble data frame with a list-col
simulating the sampling of a univoltine population with a split cohort:

- x:

  integer. x location of sampled grid cell

- y:

  integer. y location of sampled grid cell

- larvalDensity:

  integer. The density of larvae in the grid cell

- massDistribution:

  numeric. list-col of the mass of each individual in the cell

## Source

data-raw/single_cohort_sim.R
