# Sample data for community secondary production analysis

A collection of community macroinvertebrate data for estimating
secondary production

## Usage

``` r
wbtData
```

## Format

\`wbtData\` A list of two objects. 'sampleInfo' is list of 32 data
frames for each taxonomic entity and the second, 'taxaInfo' is a data
frame with 32 rows and 13 columns:

- taxonID:

  the taxonomic identifier

- repID:

  the replicate identifier

- dateID:

  the date identifier as.Date format

- lengthClass:

  the numeric (coercible) description of individual length

- n_m2:

  the count or density of individuals from a lengthClass

## Source

\<https://doi.org/10.4319/lo.2014.59.2.0507\>
