#' @description
#' This function calculates secondary production with the increment-summation method.
#' @title calc_prod_is
#' @param taxaSampleListMass description
#' @param taxaInfo data frame of taxonomic information for calculating production
#' @param bootNum integer. How many bootstrap samples should be constructed
#' @param dateDf data frame of date information sampling date such as interval length in days.
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

calc_prod_is <- function(taxaSampleListMass= NULL,
                         taxaInfo = NULL,
                         bootNum = NULL,
                         dateDf = NULL,
                         taxaSummary = TRUE,
                         wrap = FALSE,
                         massValue = 'mass',
                         abunValue = 'density',
                         dateCol = 'dateID',
                         repCol = 'repID',
                         bootList = NULL,
                         ...) {

  ## tests ##

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  # ## function prep ##
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    dateDf = dateDf,
    massValue = massValue,
    abunValue = abunValue,
    dateCol = dateCol,
    repCol = repCol,
    wrap = wrap,
    full = TRUE
  )

  # calculate the production from the full samples
  P.samp = do.call(is_prod.sample, args = funcList)
  if(P.samp$P.ann.samp == 0){
    if(taxaSummary == FALSE){
      taxaSummary <- NULL
    } else if (taxaSummary == TRUE) {
      # # create a list for output
      taxaSummary <- list(
        taxonID = speciesName,
        method = "is",
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

## perform the bootstrap procedure
  # P.boots = vector('list', length = bootNum)
  P.boots = mapply(FUN = is_prod.sample,
                   df = bootList,
                   massValue = massValue,
                   abunValue = abunValue,
                   dateCol = dateCol,
                   repCol = repCol,
                   full = FALSE,
                   MoreArgs = list(dateDf = dateDf))

# browser()
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
    taxaSummary <- NULL
  } else if (taxaSummary == TRUE) {
    # # create a list for output
    taxaSummary <- list(
      taxonID = taxaInfo$taxonID,
      method = "is",
      P.ann.samp = P.samp$P.ann.samp,
      P.uncorr.samp = NULL,
      pb = pb,
      meanN = mean(unlist(sampSummary[[paste0(abunValue,"_mean")]])),
      meanB = mean(unlist(sampSummary[["biomass_mean"]])),
      meanIndMass = P.samp$B.ann.mean / P.samp$N.ann.mean,
      datesInfo = sampSummary
    )
  }
  #   assign(taxaInfo$taxonID, list())
  # # add taxainformation
  #   assign(taxaInfo$taxonID,
  #          within(eval(as.symbol(taxaInfo$taxonID)),{
  #            taxaInfo <- taxaInfo
  #            }
  #            )
  #          )

  return(assign(speciesName, list(P.boots = P.boots,
              taxaSummary = taxaSummary)))

}
