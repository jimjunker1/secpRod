# A selection of functions to run prior to building
library(pkgbuild)
library(usethis)

# check build tool
pkgbuild::check_build_tools()

## build vignettes ----
## vignettes are included in the package

# use_vignette("simple-example", "A simple example")

## build articles ----
## articles are vignettes that are not included in the package
## but are included in the package webpage

# use_article("sampling-simulation","Simulating the sampling of populations")
# use_article("model-testing", "Testing models of secondary production: simple to advanced")


## add data ----
## to add scripts to produce package data
## these end with a call to use_data()

# use_data_raw("single_cohort_sim")

## to add additional data sets

# use_data(univoltine)
# use_data(wbtData)

# if new data objects are added,must modify R/data.R

# add author ----
usethis::use_author(
  given = "Jim",
  family = "Junker",
  role = c("aut", "cre"),
  email = "james.junker@unt.edu",
  comment = c(ORCID = "0000-0001-9713-2330")
)

# create the require documentation ----
devtools::document()

# set up the webpage ----
usethis::use_pkgdown_github_pages()

# give a package check
devtools::check()
