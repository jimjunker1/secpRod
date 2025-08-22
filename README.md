
<!-- README.md is generated from README.Rmd. Please edit that file -->

------------------------------------------------------------------------

## Overview

**secpRod** is an open source package for the analysis and calculation
of secondary production for populations and communities in R.
**secpRod** uses data from repeated sampling of population abundance and
size structure in a tall data structure and a taxon information sheet as
the base objects. In addition to estimating secondary production for
communities with multiple methods, **secpRod** also provides a tools to
perform data processing actions common in the workflow prior to
calculation, such as length to mass conversions using length~mass
relationships and visualizing population density and size structure to
aid the determination of cohort structure and assess best methods for
secondary production estimation.

## Installation

``` r
# To install the latest version Github:
# install.packages('devtools')
devtools::install_github("jimjunker1/secpRod")
```

## Status: Active Development

The `secpRod` package is currently under development for deployment in
conjunction with the ’Secondary Production and Quantitative Food Webs”
chapter in the “Methods in Stream Ecology” 5th edition. This initial
deployment will

## Usage

The basic usage of the main function, `calc_production()`, is outlined
in the [“A few simple examples”
vignette](https://jimjunker1.github.io/secpRod/articles/simple-example.html).
This article showcases the calculation of secondary production of a
[single simulated
population](https://jimjunker1.github.io/secpRod/articles/sampling-simulation.html).
For a full walkthrough, see the “A full example”

<!-- secpRod has XX functions related to data organization: -->

<!-- - `frac_merge` takes raw count data and provides areal and subsampling adjustments. This function includes QA/QC checks for un-paired fractions and returns size-class abundance standardized to area. -->

<!-- - `mass_adj` combines  -->

<!-- secpRod has XX functions related to data visualization: -->

<!-- - `len_freq` combines data returned from `frac_merge` to taxon-level histograms -->

<!-- - `allen_curve` takes  -->

<!-- `calc_production` is the main function of secpRod. This function estimates bootstrapped secondary production from multiple methods determined in the `taxaInfo` object. -->

## Contributing

This is an actively developing package and we welcome any contributions.

Please note that the ‘secpRod’ project is released with a [Contributor
Code of Conduct](.github/CODE_OF_CONDUCT.md). By contributing to this
project, you agree to abide by its terms.

If you would like to contribute please see the [contributing
document](.github/CONTRIBUTING.md).

These are heavily borrowed and adapted from the [WEEcology
lab](https://www.weecology.org/) [Portal
project](https://portal.weecology.org/) adapted from the [Contributor
Covenant](https://www.contributor-covenant.org/) and portalr package,
respectively.

## Acknowledgments

The original code base was inspired by code written by [Ben
Koch](https://scholar.google.com/citations?user=30HVNGEAAAAJ&hl=en) and
was further developed by Jim Junker with input from Wyatt F. Cross, Dan
Nelson, L. Mick Demi, and others.
