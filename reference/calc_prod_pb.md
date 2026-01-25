# calc_prod_pb

This function calculates secondary production with the
Production:Biomass method.

## Usage

``` r
calc_prod_pb(
  taxaSampleListMass = NULL,
  taxaInfo = NULL,
  bootNum = NULL,
  dateDf = NULL,
  taxaSummary = TRUE,
  wrap = TRUE,
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

  a data.frame of long format returned from
  [`convert_length_to_mass()`](https://jimjunker1.github.io/secpRod/reference/convert_length_to_mass.md)
  function

- taxaInfo:

  data frame of taxonomic information for calculating production

- bootNum:

  integer. How many bootstrap samples should be constructed

- dateDf:

  data frame of date information with external predictors for each
  month. There should be a column name identical to all variables in the
  growth equation found in taxaInfo data.frame.

- taxaSummary:

  logical. If TRUE (default) the taxaSummary will be with annual summary
  information will be returned.

- wrap:

  logical. Should the dates wrap to create a full year?

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
