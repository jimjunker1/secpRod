#' @description
#' This function calculates secondary production with the Production:Biomass method.
#' @title calc_prod_pb
#' @param taxaSampleListMass a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param taxaInfo data frame of taxonomic information for calculating production
#' @param bootNum integer. How many bootstrap samples should be constructed
#' @param dateDf data frame of date information with external predictors for each month. There should be a column name identical to all variables in the growth equation found in taxaInfo data.frame.
#' @param taxaSummary logical. If TRUE (default) the taxaSummary will be with annual summary information will be returned.
#' @param wrap logical. Should the dates wrap to create a full year?
#' @param massValue string of the column name containing the mass measurement
#' @param abunValue string of the column name containing the abundance or density measurement
#' @param dateCol string of the column name containing the date information. This can be either a recognized date object (e.g., Date, POSIX) or numeric.
#' @param repCol string of the column name containing the replicate information
#' @param bootList list. This is the bootstrapped samples passed from \code{calc_production()}
#' @param ... additional arguments to be passed to the function
#' @returns returns a list of 2 objects:
#' @returns P.boots: the boostrapped estimates of production, abundance, and biomass.
#' @returns taxaSummary: is the summary of the sample production, abundance, and biomass
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats formula
#' @export

calc_prod_pb <- function(taxaSampleListMass= NULL,
                         taxaInfo = NULL,
                         bootNum = NULL,
                         dateDf = NULL,
                         taxaSummary = TRUE,
                         wrap = TRUE,
                         massValue = 'mass',
                         abunValue = 'density',
                         dateCol = 'dateID',
                         repCol = 'repID',
                         bootList = NULL,...) {

  ## tests ##

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  taxaInfo = taxaInfo[which(taxaInfo$taxonID == speciesName),]
  # ## function prep ##
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    massValue = massValue,
    abunValue = abunValue,
    dateCol = dateCol,
    repCol = repCol,
    wrap = wrap,
    full = TRUE
  )
  # calculate the production from the full samples
  if(is.numeric(taxaInfo$pb)){
    if(length(taxaInfo$pb) == 1){
      taxaPB = unlist(taxaInfo$pb)
    } else{
      taxaPB = mean(taxaInfo$pb, na.rm = TRUE)
    }
  } else if(is.character(taxaInfo$pb)){
    taxaPB = tryCatch(
      {
        pbDist_begin = gsub("(^.*))$", "\\1", unlist(taxaInfo$pb))
        pbDist_end = gsub("^.*(\\)$)", "\\1", unlist(taxaInfo$pb))
        pbDist_n = paste0(", n = 1")
        pbDist_expr = str2expression(paste0(pbDist_begin,pbDist_n,pbDist_end))
        eval(pbDist_expr)
      },
      error = function(cond){
        message()
        return(NA_real_)
      },
      warning = function(cond){
        message(paste0("pbBoots for",speciesName," returned an error."))
      }
    )
  }

  funcList = c(funcList, list(pb = taxaPB))
  P.samp = do.call(pb_prod.sample, args = funcList)
  if(P.samp$P.ann.samp == 0){
    if(taxaSummary == FALSE){
      taxaSummary <- NULL
    } else if (taxaSummary == TRUE) {
      # # create a list for output
      taxaSummary <- list(
        taxonID = taxaInfo$taxonID,
        method = "pb",
        P.ann.samp = 0,
        pb = NA_real_,
        meanN = 0,
        meanB = 0,
        meanIndMass = 0,
        Nmean = 0,
        Nsd = 0,
        Bmean = 0,
        Bsd = 0,
        datesInfo = NULL
      )
    }
    return(assign(speciesName, list(P.boots = P.samp,
                                    taxaSummary = taxaSummary)))
  }

  # calculate the production from the full samples
  if(is.numeric(taxaInfo$pb)){
    if(length(taxaInfo$pb) == 1){
      pbBoots = rep(unlist(taxaInfo$pb), as.integer(bootNum))
    } else if(length(taxaInfo$pb) == bootNum){
      pbBoots = unlist(taxaInfo$pb)
    }} else if(is.character(taxaInfo$pb)){
      pbBoots = tryCatch(
        {
          pbDist_begin = gsub("(^.*))$", "\\1", unlist(taxaInfo$pb))
          pbDist_end = gsub("^.*(\\)$)", "\\1", unlist(taxaInfo$pb))
          pbDist_n = paste0(", n = ",as.integer(bootNum))
          pbDist_expr = str2expression(paste0(pbDist_begin,pbDist_n,pbDist_end))
          eval(pbDist_expr)
        },
        error = function(cond){
          message()
          return(NA_real_)
        },
        warning = function(cond){
          message(paste0("pbBoots for",speciesName," returned an error."))
        }
      )
  }

  P.boots = mapply(FUN = pb_prod.sample,
                   df = bootList,
                   massValue = massValue,
                   abunValue = abunValue,
                   dateCol = dateCol,
                   repCol = repCol,
                   pb = pbBoots,
                   full = FALSE)
  #### create SAMPLE information to export as summary ####
  sampSummary = create_sample_summary(df = taxaSampleListMass,
                                      wrap = wrap,
                                      abunValue = abunValue,
                                      massValue = massValue,
                                      dateCol = dateCol,
                                      repCol = repCol,
                                      ...)

  #estimate the sample PB
  pb = P.samp$P.ann.samp/P.samp$B.ann.mean
  if(taxaSummary == FALSE){

  } else if (taxaSummary == TRUE) {
    # # create a list for output
    taxaSummary <- list(
      taxonID = taxaInfo$taxonID,
      method = "pb",
      P.ann.samp = P.samp$P.ann.samp,
      pb = pb,
      meanN = P.samp$N.ann.mean,
      meanB = P.samp$B.ann.mean,
      meanIndMass = P.samp$B.ann.mean /P.samp$N.ann.mean,
      datesInfo = sampSummary
    )
  }

  return(assign(speciesName, list(P.boots = P.boots,
              taxaSummary = taxaSummary)))

}
