#' @name calc_production
#' @title calc_production
#' @description This is the main function of the secpRod package. It will calculate secondary production for all groups based on the methods described in the taxa information object. Depending on input values varying summaries are returned.
#' @param taxaSampleListMass data.frame of sample length-masses and abundances
#' @param infoCols integer vector; Any columns in the sizeInfo object that are not the taxonomic ID, sampling metadata, or size class columns. ** This is soft-deprecated and is not currently used. It may go away soon **
#' @param taxaInfo dataframe (or coercible); The taxonomic information for calculating secondary production. This must include a taxonomic ID column with the same name as that of \code{taxaSampleListMass}
#' @param bootNum integer. The number of bootstrapped samples to create
#' @param taxaSummary logical. If TRUE (default) the taxaSummary will be with annual summary information will be returned.
#' @param lengthValue string of the column name containing length class measurements. If NULL (default)
#' @param massValue string of the column name containing the mass measurement
#' @param abunValue string of the column name containing the abundance or density measurement
#' @param dateCol string of the column name containing the date information. This can be either a recognized date object (e.g., Date, POSIX)
#' @param repCol string of the column name containing the replicate information
#' @param envData data.frame of additional environmental data for each sampling event.  If present, this is combined with the dateDf information and passed to the IGR function for growth rate estimations.The environmental variables passed should align with the start of each sampling interval.
#' @param ... additional arguments passed to function
#'@export

calc_production = function(taxaSampleListMass = NULL,
                           # infoCols = NULL,
                           taxaInfo = NULL,
                           bootNum = 1e2,
                           taxaSummary = TRUE,
                           lengthValue = NULL,
                           massValue = 'mass',
                           abunValue = 'density',
                           dateCol = 'dateID',
                           repCol = 'repID',
                           envData = NULL,
                           ...){

  ### tests ###
  # taxaInfo is currently only allowed to have one species.
  # This function is currently only set up for a single species. Future updates will allow full community
  # estimates to be done with a single call.
  if(length(unlist(taxaInfo$taxonID)) > 1) stop("Error: More than one species' taxaInfo passed to function. Only single species are allowed within each call currently.")
  # are all the input values in the taxaSampleList?
  colVec = c(lengthValue, massValue, abunValue, dateCol, repCol)
  if(!all(colVec %in% names(taxaSampleListMass))){
    stop(paste0("Error: Not all named column inputs found in the taxaSampleListMass. ", paste(colVec[colVec %ni% names(taxaSampleListMass)], collapse = ",")," were not present."))
  }
  # are the methods all recognizable?
  if(!all(unlist(taxaInfo$method) %in% c('is','sf','pb','igr'))){
    badMethod = unique(unlist(taxaInfo$method)[which(unlist(taxaInfo$method) %ni% c('is','sf','pb','igr'))])
    stop(paste0("Error: ",badMethod," is not a recognized method. Available values are 'is','sf','pb','igr'. See documentation for more information."))
  }
  # is the date object a POSIX? convert to Date for merging.
  if(inherits(taxaSampleListMass[[dateCol]], c("POSIXt"))){
    taxaSampleListMass[[dateCol]] <- as.Date(taxaSampleListMass[[dateCol]])
  }
  ### end tests ###
  # prep size-abundance bootstraps. these are fixed across methods if multiple methods passed for comparisons
  bootList = vector('list', length = bootNum)
  bootList = prep_boots(df = taxaSampleListMass,
                        bootNum = bootNum,
                        dateCol = dateCol,
                        repCol = repCol)
  ## function prep ##
  speciesName = unique(taxaSampleListMass$taxonID)

  ### make a list of key variables to pass
  massValue = massValue
  abunValue = abunValue

  # need to make this into method specific if multiple methods are passed?
  wrap = unlist(taxaInfo$wrap)
  # build the dateDf object by combining wrap_dates with envData
  # suppressWarnings of coercing double to logical
  dateDf = suppressWarnings(wrap_dates(df = taxaSampleListMass, envData = envData, dateCol = dateCol, wrapDate = wrap))

  # build the sizesDf object depending on if lengthValue is present or not
  if(is.null(lengthValue)){
    sizesDf = unique(taxaSampleListMass[, c(massValue)])
  } else{
    sizesDf = unique(taxaSampleListMass[, c(lengthValue, massValue)])
  }

  # create the funcList object to pass to all production methods
  funcList <- list(
    taxaSampleListMass = taxaSampleListMass,
    # infoCols = infoCols,
    taxaInfo = taxaInfo,
    lengthValue = lengthValue,
    massValue = massValue,
    abunValue = abunValue,
    dateCol = dateCol,
    repCol = repCol,
    dateDf = dateDf,
    # dataframe of sizes and masses
    sizesDf = sizesDf,
    bootNum = bootNum,
    taxaSummary = taxaSummary,
    wrap = wrap,
    bootList = bootList
  )
### determine method specific funcList variables ###
# Are there multiple methods in nested list-col? Unnest them and add to funcList
  if('list' %in% class(taxaInfo$method)){
    funcList[['method']] = list(unlist(taxaInfo$method))
  } else{
    funcList[['method']] = unlist(taxaInfo$method)
  }

# are there notes in the taxaInfo object?
  if(!is.null(taxaInfo$notes)){
    notes = taxaInfo$notes
  } else notes = NULL


## calculated production based on methods

### increment summation
  if('is' %in% funcList$method){
    is_prod = do.call(calc_prod_is, args = funcList)
    return(is_prod)
  }
### size frequency
  if('sf' %in% funcList$method){
  sf_prod = do.call(calc_prod_sf, args = funcList)
  return(sf_prod)
  }
  # } else if(taxaInfo$method == "igr"){
### instantaneous growth
  if('igr' %in% funcList$method){
  igr_prod = do.call(calc_prod_igr, args = funcList)
  return(igr_prod)
  }
### production:biomass method
  if('pb' %in% taxaInfo$method){
  pb_prod = do.call(calc_prod_pb, args =funcList)
  return(pb_prod)
  }
### removal summation
  # rs_prod = calc_prod_rs()


}
