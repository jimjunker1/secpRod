# calc_prod_igr

This function calculates production using the instantaneous growth rate
method.

## Usage

``` r
calc_prod_igr(
  taxaSampleListMass = NULL,
  taxaInfo = NULL,
  bootNum = NULL,
  dateDf = NULL,
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

  The taxa info data.frame

- bootNum:

  integer. The number of bootstrap samples to produce.

- dateDf:

  data frame of date information with external predictors for each
  month. There should be a column name identical to all variables in the
  growth equation found in taxaInfo data.frame.

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

## Details

Additional details...
