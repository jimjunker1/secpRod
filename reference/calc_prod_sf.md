# calc_prod_sf

This function calculates secondary production with the size-frequency
method.

## Usage

``` r
calc_prod_sf(
  taxaSampleListMass = NULL,
  taxaInfo = NULL,
  bootNum = NULL,
  taxaSummary = TRUE,
  wrap = FALSE,
  lengthValue = NULL,
  massValue = "mass",
  abunValue = "density",
  dateCol = "dateID",
  repCol = "repID",
  bootList = NULL,
  ...
)
```

## Arguments

- taxaSampleListMass:

  description

- taxaInfo:

  data frame of taxonomic information for calculating production

- bootNum:

  integer. How many bootstrap samples should be constructed

- taxaSummary:

  logical. If TRUE (default) the taxaSummary will be with annual summary
  information will be returned.

- wrap:

  logical. Should the dates wrap to create a full year?

- lengthValue:

  string of the column name containing the length class measurements

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

- bootList:

  list. This is the bootstrapped samples passed from
  [`calc_production()`](https://jimjunker1.github.io/secpRod/reference/calc_production.md)

- ...:

  additional arguments to be passed to the function

## Value

returns a list of 2 objects:

P.boots: the boostrapped estimates of production, abundance, and
biomass.

taxaSummary: is the summary of the sample production, abundance, and
biomass
