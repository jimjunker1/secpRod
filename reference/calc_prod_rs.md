# calc_prod_rs

This function calculates secondary production with the removal-summation
method.

## Usage

``` r
calc_prod_rs(
  taxaSampleListMass = NULL,
  taxaInfo = NULL,
  bootNum = NULL,
  dateDf = NULL,
  taxaSummary = "full",
  wrap = TRUE,
  massValue = NULL,
  massLabel = NULL,
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

- dateDf:

  data frame of date information with external predictors for each
  month. There should be a column name identical to all variables in the
  growth equation found in taxaInfo data.frame.

- taxaSummary:

  string of 'short', 'full', or 'none'. What type of summary information
  should be returned.

- wrap:

  logical. Should the dates wrap to create a full year?

- massValue:

  string. What is the mass value and units of the production

- massLabel:

  string. What label should the output units be. It is possible this
  will default to 'massValue' in the future.

- bootList:

  list. This is the bootstrapped samples passed from
  \`calc_production()\`

- ...:

  additional arguments to be passed to the function

## Value

returns a list of 2 objects:

P.boots: the boostrapped estimates of production, abundance, and
biomass.

taxaSummary: is the summary of the sample production, abundance, and
biomass
