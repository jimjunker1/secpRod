# plot_cohorts

This function plots size-frequency data of a single taxon over time

## Usage

``` r
plot_cohorts(
  taxaSampleListMass = NULL,
  param = c("length", "mass"),
  lengthValue = NULL,
  massValue = "mass",
  abunValue = NULL,
  dateCol = "dateID",
  repCol = "repID",
  ...
)
```

## Arguments

- taxaSampleListMass:

  a data.frame of long format returned from
  [`convert_length_to_mass()`](https://jimjunker1.github.io/secpRod/reference/convert_length_to_mass.md)
  function

- param:

  character. A string of 'length' or 'mass' that describes what
  measurement should be plotted

- lengthValue:

  a character string idenfitying the column name with the length values
  / length classes

- massValue:

  character. The column name of the mass measurement (e.g., mass,
  afdmMg, etc.)

- abunValue:

  a character string identifying the column name with the abundance
  values.

- dateCol:

  a character string identifying the column name with the date
  information

- repCol:

  a character string identifying the column name with the replicate
  information

- ...:

  additional arguments passed to function

## Value

returns a histogram of the plot of the relative frequency of size or
mass classes for a single taxon for all sampling dates
