#' @description
#' This function calculates secondary production with the increment-summation method.
#' @title calc_prod_is
#' @param taxaSampleListMass description
#' @param taxaInfo data frame of taxonomic information for calculating production
#' @param bootNum integer. How many bootstrap samples should be constructed
#' @param dateDf data frame of date information with external predictors for each month. There should be a column name identical to all variables in the growth equation found in taxaInfo data.frame.
#' @param taxaSummary string of 'short', 'full', or 'none'. What type of summary information should be returned.
#' @param wrap logical. Should the dates wrap to create a full year?
#' @param massValue string. What is the mass value and units of the production
#' @param massLabel string. What label should the output units be. It is possible this will default to 'massValue' in the future.
#' @param bootList list. This is the bootstrapped samples passed from `calc_production()`
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
                         taxaSummary = 'full',
                         wrap = FALSE,
                         massValue = 'afdm_mg',
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
    if(taxaSummary == "none"){
      taxaSummary <- NULL
    } else if (taxaSummary == "full") {
      # # create a list for output
      taxaSummary <- list(
        summaryType = "full",
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
    } else if(taxaSummary == "short"){
      taxaSummary <- list(
        summaryType = "short",
        taxonID = taxaInfo$taxonID,
        method = "is",
        P.ann.samp = 0,
        pb = NA_real_,
        meanN = 0,
        meanB = 0,
        meanIndMass = 0,
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
  if(taxaSummary == "none"){

  } else if (taxaSummary == "full") {
    # # create a list for output
    taxaSummary <- list(
      summaryType = "full",
      taxonID = taxaInfo$taxonID,
      method = "is",
      P.ann.samp = P.samp$P.ann.samp,
      P.uncorr.samp = NULL,
      # cpi = taxaCPI,
      pb = pb,
      meanN = mean(unlist(sampSummary[[paste0(abunValue,"_mean")]])),
      meanB = mean(unlist(sampSummary[["biomass_mean"]])),
      meanIndMass = P.samp$B.ann.mean / P.samp$N.ann.mean,
      # Nmean = NmeanTab,
      # Nsd = NsdTab,
      # Bmean = BmeanTab,
      # Bsd = BsdTab,
      datesInfo = sampSummary
    )
  } else if(taxaSummary == "short"){
    taxaSummary <- list(
      summaryType = "short",
      taxonID = taxaInfo$taxonID,
      method = "is",
      P.ann.samp = P.samp$P.ann.samp,
      # cpi = taxaCPI,
      pb = pb,
      meanN = mean(unlist(sampSummary$n_m2_mean)),
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
