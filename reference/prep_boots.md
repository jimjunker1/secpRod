# prep_boots

This function prepares bootstrap samples.

## Usage

``` r
prep_boots(df = NULL, bootNum = bootNum, dateCol = "dateID", repCol = "repID")
```

## Arguments

- df:

  data.frame. A dataframe of the species size, mass, and frequency data.

- bootNum:

  integer. The number of bootstrapped data sets that should be created.

- dateCol:

  a string of the column name containing date information

- repCol:

  a string of the column name containing the replicate information
