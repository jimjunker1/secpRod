# rs_prod.sample

This function calculates taxa production based on the removal-summation
method

## Usage

``` r
rs_prod.sample(
  df = NULL,
  sizesDf = NULL,
  massValue = NULL,
  massLabel = NULL,
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

  a data.frame of the

- massValue:

  character string identifying the column name of the mass value

- massLabel:

  character string identifying the column name of the mass value

- full:

  logical. should the full summary be returned with mean and sd

- ...:

  additional arguments passed to function

## Value

list object with taxa summary of the sampled data
