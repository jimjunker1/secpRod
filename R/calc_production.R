#' @name calc_production
#' @title calc_production
#' @description This is the main function of the secpRod package. It will calculate secondary production for all groups based on the methods described in the taxa information object. Depending on input values varying summaries are returned.
#' @param taxaSampleListMass data.frame of sample length-masses and abundances
#' @param infoCols integer vector; Any columns in the sizeInfo object that are not the taxonomic ID, sampling metadata, or size class columns
#' @param taxaInfo dataframe (or coercible); The taxonomic information for calculating secondary production. This must include a taxonomic ID column with the same name as that of \code{taxaSampleListMass}
#' @param bootNum integer. The number of bootsrapt samples to create
#' @param wrap logical should an extra date be added to make a full calendar year?
#' @param taxaSummary string of \code{'short'}, \code{'full'}, or \code{'none'} to distinguish the information returned
#' @param ... additional arguments passed to function
#'@export

calc_production = function(taxaSampleListMass = NULL,
                           infoCols = NULL,
                           taxaInfo = NULL,
                           bootNum = 1e2,
                           wrap = 1L,
                           taxaSummary = 'full',
                           ...){
  ### tests ###

  ### end tests ###
  ## function prep ##
  ##
  ### make a list of key variables to pass
  massValue = rev(names(taxaSampleListMass))[1]
  massLabel = paste0(massValue, "_m2")

  funcList <- list(
    taxaSampleListMass = taxaSampleListMass,
    infoCols = infoCols,
    taxaInfo = taxaInfo,
    massValue = massValue,
    massLabel = massLabel,
    dateDf = wrap_dates(df = taxaSampleListMass, wrapDate = wrap),
    # dataframe of sizes and masses
    sizesDf = unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])]),
    bootNum = bootNum,
    taxaSummary = taxaSummary,
    wrap = TRUE
  )

## calculated production based on methods
### size frequency
#   debugonce(sf_prod.sample)
#   debugonce(calc_prod_sf)
  sf_prod = do.call(calc_prod_sf, args = funcList)
### instantaneous growth
  # igr_prod = do.call(calc_prod_igr, args = funcList)
### removal summation
  # rs_prod = calc_prod_rs()
### increment summation
  # is_prod = calc_prod_is()

  return(sf_prod)

}
