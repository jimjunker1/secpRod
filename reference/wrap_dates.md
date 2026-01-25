# wrap_dates

This function creates an additional date in the data set to make a full
annual cycle

## Usage

``` r
wrap_dates(df = NULL, envData = NULL, dateCol = NULL, wrapDate = TRUE)
```

## Arguments

- df:

  a data.frame of long format returned from
  [`convert_length_to_mass()`](https://jimjunker1.github.io/secpRod/reference/convert_length_to_mass.md)
  function

- envData:

  (optional) a data.frame of the environmental data. The number of rows
  should be the same as the number of sampling dates.

- dateCol:

  character string identifying the column name of the sampling dates

- wrapDate:

  logical. should we wrap the dates to make full year?

## Value

data.frame with sampling dates and intervals, or number of days between
