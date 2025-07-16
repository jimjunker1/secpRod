#' Sample data for community secondary production analysis
#'
#' A collection of community macroinvertebrate data for estimating secondary production
#'
#' @format `wbtData`
#' A list of two objects. 'sampleInfo' is list of 32 data frames for each taxonomic entity and the second, 'taxaInfo' is a data frame with 32 rows and 13 columns:
#' \describe{
#'   \item{taxonID}{the taxonomic identifier}
#'   \item{repID}{the replicate identifier}
#'   \item{dateID}{the date identifier as.Date format}
#'   \item{lengthClass}{the numeric (coercible) description of individual length}
#'   \item{n_m2}{the count or density of individuals from a lengthClass}
#' }
#' @source <https://doi.org/10.4319/lo.2014.59.2.0507>
"wbtData"


#' Sample data for community secondary production analysis
#'
#' A collection of community macroinvertebrate data for estimating secondary production
#'
#' @format `univoltine`
#' A single data.frame of the sampleInfo data as sampled for a single simulated taxon with columns:
#' \describe{
#'   \item{taxonID}{the taxonomic identifier}
#'   \item{repID}{the replicate identifier}
#'   \item{dateID}{the date identifier as.Date format}
#'   \item{lengthClass}{the numeric (coercible) description of individual length}
#'   \item{n_m2}{the count or density of individuals from a lengthClass}
#' }
"univoltine"

#' Sample data for community secondary production analysis
#'
#' A collection of community macroinvertebrate data for estimating secondary production
#'
#' @format `singleCohortSim`
#' A single tibble data frame with a list-col simulating the sampling of a single cohort population:
#' \describe{
#'   \item{x}{integer. x location of sampled grid cell}
#'   \item{y}{integer. y location of sampled grid cell}
#'   \item{larvalDensity}{integer. The density of larvae in the grid cell}
#'   \item{massDistribution}{numeric. list-col of the mass of each individual in the cell}
#' }
#' @source data-raw/single_cohort_sim.R
"singleCohortSim"


#' Sample data for community secondary production analysis
#'
#' A collection of community macroinvertebrate data for estimating secondary production
#'
#' @format `splitCohortSim`
#' A single tibble data frame with a list-col simulating the sampling of a univoltine population with a split cohort:
#' \describe{
#'   \item{x}{integer. x location of sampled grid cell}
#'   \item{y}{integer. y location of sampled grid cell}
#'   \item{larvalDensity}{integer. The density of larvae in the grid cell}
#'   \item{massDistribution}{numeric. list-col of the mass of each individual in the cell}
#' }
#' @source data-raw/single_cohort_sim.R
"splitCohortSim"

#' Sample data for community secondary production analysis
#'
#' A collection of community macroinvertebrate data for estimating secondary production
#'
#' @format `overlapCohortSim`
#' A single tibble data frame with a list-col simulating the sampling of a univoltine population with a split cohort:
#' \describe{
#'   \item{x}{integer. x location of sampled grid cell}
#'   \item{y}{integer. y location of sampled grid cell}
#'   \item{larvalDensity}{integer. The density of larvae in the grid cell}
#'   \item{massDistribution}{numeric. list-col of the mass of each individual in the cell}
#' }
#' @source data-raw/single_cohort_sim.R
"overlapCohortSim"
