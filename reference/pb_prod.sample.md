# pb_prod.sample

This function calculates sample production based on the
production:biomass method

## Usage

``` r
pb_prod.sample(
  df = NULL,
  massValue = "mass",
  abunValue = "density",
  dateCol = "dateID",
  repCol = "repID",
  pb = NULL,
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

- massValue:

  character string identifying the column name of the mass value

- abunValue:

  character string identifying the column name of the abundance value

- dateCol:

  character string identifying the column name of the date column

- repCol:

  character string identifying the column name of the replicate column

- pb:

  numeric. The production:biomass ratio.

- wrap:

  logical. should the calculations be wrapped by adding an additional
  date to make a full year?

- full:

  logical. should the full summary be returned with mean and sd

- ...:

  additional arguments passed to function

## Value

list object with annual production, mean biomass, and mean abundance of
the sample
