# create_sample_summary

A function to convert species-specific lengths to mass based on a
user-provided length-mass equation form and variable values

## Usage

``` r
create_sample_summary(
  df = NULL,
  wrap = FALSE,
  abunValue = NULL,
  massValue = NULL,
  dateCol = NULL,
  repCol = NULL,
  ...
)
```

## Arguments

- df:

  data.frame in long format for a single taxa. The data.frame should
  contain a species identifier column \`taxonID\` and a column of length
  bin categories \`lengthClass\`. \`lengthClass\` values must be numeric
  or coercible.

- wrap:

  logical (or coercible) indicating whether the production estimate
  should wrap the first and last dates to create a full annual cycle. If
  TRUE, this will create an additional sampling interval using the mean
  densities and masses to create a full annual data set.

- abunValue:

  string of the column name containing the abundance or density
  measurement

- massValue:

  string of the column name containing the mass measurement

- dateCol:

  string of the column name containing the date information. This can be
  either a recognized date object (e.g., Date, POSIX)

- repCol:

  string of the column name containing the replicate information

- ...:

  additional arguments passed to function

## Value

a summary data.frame of the sample summary for N: sample size, density:
mean and sd, biomass: mean and sd
