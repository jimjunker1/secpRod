# A selection of functions to run prior to building
library(pkgbuild)
library(usethis)

# check build tool
pkgbuild::check_build_tools()

# build vignettes ----
use_article("simple-example", "A simple example")



# add data ----
use_data() # e.g. use_data(univoltine)

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
