#' @name calc_production
#' @title calc_production
#' @description This is the main function of the secpRod package. It will calculate secondary production for all groups based on the methods described in the taxa information object. Depending on input values varying summaries are returned.
#' @param taxaSampleListMass data.frame of sample length-masses and abundances
#' @param infoCols integer vector; Any columns in the sizeInfo object that are not the taxonomic ID, sampling metadata, or size class columns
#' @param taxaInfo dataframe (or coercible); The taxonomic information for calculating secondary production. This must include a taxonomic ID column with the same name as that of \code{taxaSampleListMass}
#' @param bootNum integer. The number of bootstrapped samples to create
#' @param taxaSummary string of \code{'short'}, \code{'full'}, or \code{'none'} to distinguish the information returned
#' @param lengthValue string of the column name containing length class measurements. If NULL (default)
#' @param massValue string of the column name containing the mass measurement
#' @param abunValue string of the column name containing the abundance or density measurement
#' @param ... additional arguments passed to function
#'@export

calc_production = function(taxaSampleListMass = NULL,
                           infoCols = NULL,
                           taxaInfo = NULL,
                           bootNum = 1e2,
                           taxaSummary = 'full',
                           lengthValue = NULL,
                           massValue = 'afdm_mg',
                           abunValue = 'density',
                           ...){

  ### tests ###
  ## taxaInfo is currently only allowed to have one species.
  ## This function is currently only set up for a single species. Future updates will allow full community
  ## estimates to be done with a single call.
  if(length(droplevels(unlist(taxaInfo$taxonID))) > 1) stop("Error: More than one species' taxaInfo passed to function. Only single species are allowed within each call currently.")

  ### end tests ###
  # prep size-abundance boots
  bootList = prep_boots(df = taxaSampleListMass,
                        bootNum = bootNum)
  ## function prep ##
  speciesName = unique(taxaSampleListMass$taxonID)

  ### make a list of key variables to pass
  massValue = massValue
  abunValue = abunValue
  # massLabel = paste0(massValue, "_m2")
  wrap = unlist(taxaInfo$wrap)

  funcList <- list(
    taxaSampleListMass = taxaSampleListMass,
    infoCols = infoCols,
    taxaInfo = taxaInfo,
    massValue = massValue,
    # massLabel = massLabel,
    dateDf = wrap_dates(df = taxaSampleListMass, wrapDate = wrap),
    # dataframe of sizes and masses
    sizesDf = unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])]),
    bootNum = bootNum,
    taxaSummary = taxaSummary,
    wrap = wrap,
    bootList = bootList
  )

## calculated production based on methods

### size frequency
  if('sf' %in% taxaInfo$method){
#   debugonce(sf_prod.sample)
#   debugonce(calc_prod_sf)
  sf_prod = do.call(calc_prod_sf, args = funcList)
  return(sf_prod)
  # } else if(taxaInfo$method == "igr"){
### instantaneous growth
  # igr_prod = do.call(calc_prod_igr, args = funcList)
  } else if('pb' %in% taxaInfo$method){
### pb
  pb_prod = do.call(calc_prod_pb, args =funcList)
  return(pb_prod)
  }
### removal summation
  # rs_prod = calc_prod_rs()
### increment summation
  # is_prod = calc_prod_is()

  # return(sf_prod)
  # return(pb_prod)

}
