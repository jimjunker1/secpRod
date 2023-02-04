#' @name calc_production
#' @description This is the main function of the secpRod package. It will calculate secondary production for all groups based on the methods described in the taxa information object. Depending on input values varying summaries
#' @param sizeInfo dataframe (or coercible); The size-abundance data for each species and sampling date
#' @param infoCols integer vector; Any columns in the sizeInfo object that are not the taxonomic ID, sampling data, or size class columns
#' @param taxaInfo dataframe (or coercible); The taxonomic information for calculating secondary production. This must include a taxonomic ID column with the same name as that of `sizeInfo`, other information
#'

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

  debugonce(calc_prod_sf)
  do.call(calc_prod_sf, args = funcList)
  # calc_prod_sf(df = taxaSampleListMass,  )
  # calc_prod_rs()
  # calc_prod_is()
  # calc_prod_igr()

}
