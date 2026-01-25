# calc_production

This is the main function of the secpRod package. It will calculate
secondary production for all groups based on the methods described in
the taxa information object. Depending on input values varying summaries
are returned.

## Usage

``` r
calc_production(
  taxaSampleListMass = NULL,
  taxaInfo = NULL,
  bootNum = 100,
  taxaSummary = TRUE,
  lengthValue = NULL,
  massValue = "mass",
  abunValue = "density",
  dateCol = "dateID",
  repCol = "repID",
  envData = NULL,
  ...
)
```

## Arguments

- taxaSampleListMass:

  data.frame of sample length-masses and abundances

- taxaInfo:

  dataframe (or coercible); The taxonomic information for calculating
  secondary production. This must include a taxonomic ID column with the
  same name as that of `taxaSampleListMass`

- bootNum:

  integer. The number of bootstrapped samples to create

- taxaSummary:

  logical. If TRUE (default) the taxaSummary will be with annual summary
  information will be returned.

- lengthValue:

  string of the column name containing length class measurements. If
  NULL (default)

- massValue:

  string of the column name containing the mass measurement

- abunValue:

  string of the column name containing the abundance or density
  measurement

- dateCol:

  string of the column name containing the date information. This can be
  either a recognized date object (e.g., Date, POSIX)

- repCol:

  string of the column name containing the replicate information

- envData:

  data.frame of additional environmental data for each sampling event.
  If present, this is combined with the dateDf information and passed to
  the IGR function for growth rate estimations.The environmental
  variables passed should align with the start of each sampling
  interval.

- ...:

  additional arguments passed to function
