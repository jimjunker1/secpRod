
<!-- README.md is generated from README.Rmd. Please edit that file -->

-----

## Overview

**secpRod** is an open source package for the analysis and calculation
of secondary production for populations and communities in R.
**secpRod** uses data from repeated sampling of population abundance and
size structure in a tall data structure and a taxon information sheet as
the base objects. In addition to estimating secondary production for
communities with multiple methods, **secpRod** also allows the user to
visualize population data to assess best methods for secondary
production estimation.

## Installation

``` r
# To install the latest version Github:
# install.packges('devtools')
devtools::install_github("jimjunker1/secpRod")
```

## Functions

secpRod has XX functions related to data organization:

  - `frac_merge` takes raw count data and provides areal and subsampling
    adjustments. This function includes QA/QC checks for un-paired
    fractions and returns size-class abundance standardized to area.
  - `mass_adj` combines

secpRod has XX functions related to data visualization:

  - `len_freq` combines data returned from `frac_merge` to taxon-level
    histograms
  - `allen_curve` takes

-----

Please note that the ‘secpRod’ project is released with a [Contributor
Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

-----
