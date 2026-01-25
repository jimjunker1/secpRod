# is_prod.sample

This function calculates taxa production based on the
increment-summation method

## Usage

``` r
is_prod.sample(
  df = NULL,
  dateDf = dateDf,
  massValue = "mass",
  abunValue = "density",
  dateCol = "dateID",
  repCol = "repID",
  wrap = FALSE,
  full = TRUE,
  ...
)
```

## Arguments

- df:

  a data.frame of long format returned from
  [`convert_length_to_mass()`](https://jimjunker1.github.io/secpRod/reference/convert_length_to_mass.md)
  function

- dateDf:

  data frame of date information sampling date such as interval length
  in days.

- massValue:

  string of the column name containing the mass measurement

- abunValue:

  string of the column name containing the abundance or density
  measurement

- dateCol:

  string of the column name containing the date information. This can be
  either a recognized date object (e.g., Date, POSIX) or numeric.

- repCol:

  string of the column name containing the replicate information

- wrap:

  logical. should the calculations be wrapped by adding an additional
  date to make a full year?

- full:

  logical. should the summary be returned with mean and sd

- ...:

  additional arguments passed to function

## Value

list object with annual production, mean biomass, and mean abundance of
the sample
