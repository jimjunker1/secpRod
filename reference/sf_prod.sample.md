# sf_prod.sample

This function calculates taxa production based on the size-frequency
method

## Usage

``` r
sf_prod.sample(
  df = NULL,
  sizesDf = NULL,
  lengthValue = NULL,
  massValue = "mass",
  abunValue = "density",
  dateCol = "dateID",
  repCol = "repID",
  wrap = FALSE,
  cpi = NULL,
  full = TRUE,
  ...
)
```

## Arguments

- df:

  a data.frame of long format returned from
  [`convert_length_to_mass()`](https://jimjunker1.github.io/secpRod/reference/convert_length_to_mass.md)
  function

- sizesDf:

  a data.frame of the size class including lengthClass, massClass,
  bin_min, bin_max, midpoint

- lengthValue:

  string of the column name containing the length class measurements

- massValue:

  character string identifying the column name of the mass value

- abunValue:

  character string identifying the column name of the density value

- dateCol:

  character string identifying the column name of the sample date
  information. This is distinguished from \*Value parameters in that
  Values may be used to maintain units provenance in the future. This
  may also change.

- repCol:

  character string identifying the column name of the replicate
  information.

- wrap:

  logical. should the calculations be wrapped by adding an additional
  date to make a full year?

- cpi:

  integer. The cohort production interval.

- full:

  logical. should the full summary be returned with mean and sd

- ...:

  additional arguments passed to function

## Value

list object with annual production, mean biomass, and mean abundance of
the sample
