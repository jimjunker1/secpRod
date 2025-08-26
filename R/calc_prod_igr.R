#' @description
#' This function calculates production using the instantaneous growth rate method.
#' @details
#' Additional details...
#' @title calc_prod_igr
#' @param taxaSampleListMass description
#' @param taxaInfo The taxa info data.frame
#' @param bootNum integer. The number of bootstrap samples to produce.
#' @param dateDf data frame of date information with external predictors for each month. There should be a column name identical to all variables in the growth equation found in taxaInfo data.frame.
#' @param taxaSummary logical. If TRUE (default) the taxaSummary will be with annual summary information will be returned.
#' @param wrap logical. Should the dates wrap to create a full year?
#' @param lengthValue string of the column name containing the length class measurements
#' @param massValue string of the column name containing the mass measurement
#' @param abunValue string of the column name containing the abundance or density measurement
#' @param dateCol string of the column name containing the date information. This can be either a recognized date object (e.g., Date, POSIX) or numeric.
#' @param repCol string of the column name containing the replicate information
#' @param bootList list. This is the bootstrapped samples passed from \code{calc_production()}
#' @param ... additional arguments to be passed to the function
#' @returns returns a list of 2 objects:
#' @returns P.boots: the boostrapped estimates of production, abundance, and biomass.
#' @returns taxaSummary: is the summary of the sample production, abundance, and biomass
#' @importFrom stats formula aggregate na.omit
#' @importFrom formula.tools lhs rhs
#' @importFrom rlang := sym
#' @importFrom tidyr unnest
#' @export

calc_prod_igr <- function(taxaSampleListMass= NULL,
                         taxaInfo = NULL,
                         bootNum = NULL,
                         dateDf = NULL,
                         taxaSummary = TRUE,
                         wrap = FALSE,
                         lengthValue = NULL,
                         massValue = 'mass',
                         abunValue = 'density',
                         dateCol = 'dateID',
                         repCol = 'repID',
                         bootList = NULL,
                         ...) {

  ## tests ##
  speciesName = unique(taxaSampleListMass$taxonID)

  # make sure all the variables in the growth formula are in dateDf or taxaSampleListMass
  ## remove white space to each parsing
  growthFormula = stats::formula(gsub(" ","", taxaInfo$growthForm))

  # get the left hand side and confirm "g_d" is present
  growthUnits <- formula.tools::lhs(growthFormula)
  if(!any(grepl('g_d', growthUnits))) stop("Error: the left hand side of `growthForm` must contain 'g_d'.")
  growthLHS <- deparse(formula.tools::lhs(growthFormula))
  # get the right hand side and confirm all variables in dateDf or taxaSampleListMass
  growthRHS <- formula.tools::rhs(growthFormula)
  # remove code calls, e.g. +, -, ^ , log(), etc.
  charRHS <- deparse(growthRHS)
  # split on punctuation to extract just variable names
  allGrowthVars <- sapply(unlist(strsplit(charRHS, "[[:punct:]]")), trimws)
  # remove any blank spaces from the vector
  allGrowthVars <- allGrowthVars[grepl("\\w", allGrowthVars)]
  # remove any numeric
  allGrowthCharVars = na.omit(allGrowthVars[(suppressWarnings(is.na(as.numeric(allGrowthVars))))])
  # remove "log" if present
  allGrowthCharVars = allGrowthCharVars[!grepl("log|exp",allGrowthCharVars)]
  # check that all growth variables are present in the data
  allDataVars <- c(names(taxaSampleListMass),names(dateDf))
  anyMissing <- allGrowthCharVars %ni% allDataVars
  if(sum(anyMissing) > 0) stop(paste0("Error: ",paste(allGrowthCharVars[anyMissing], collapse = ",")," are missing from sample info and environmental data."))
  ## end tests ##
  ## function prep ##
  ### create a data.frame of growth rates by size and date
  ### remove all size classes with density of 0
  taxaSampleListMassNOzeros <- dplyr::filter(taxaSampleListMass, !!rlang::sym(abunValue) > 0)
  sizeValue <- allGrowthCharVars[allGrowthCharVars %in% c(lengthValue, massValue)]
  sizeDateAggForm <- paste0(sizeValue,"~",dateCol)
  sizesDf <- tidyr::unnest(stats::aggregate(stats::formula(sizeDateAggForm), data = taxaSampleListMassNOzeros, FUN = unique), cols = !!rlang::sym(sizeValue))
  ### merge with dateDf to get date, interval length and envData
  growthFormula <- stats::as.formula(paste(growthLHS,"~",charRHS, collapse = " "))
  growthDf <- merge(dateDf, sizesDf, by = dateCol)
  growthDf <- dplyr::mutate(growthDf, !!growthLHS := !!growthRHS)
  if(!is.null(taxaInfo$min.growth)){
    growthDf[growthDf$g_d < 0, "g_d"] <- taxaInfo$min.growth
  }
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    dateDf = dateDf,
    growthDf = growthDf,
    lengthValue = lengthValue,
    massValue = massValue,
    abunValue = abunValue,
    dateCol = dateCol,
    repCol = repCol,
    wrap = wrap
  )

  # calculate the production from the observed samples
  P.samp = do.call(igr_prod.sample, args = funcList)
  if(P.samp$P.ann.samp == 0){
    if(taxaSummary == FALSE){
      taxaSummary <- NULL
    } else if (taxaSummary == TRUE) {
      # # create a list for output
      taxaSummary <- list(
        taxonID = taxaInfo$taxonID,
        method = "igr",
        P.ann.samp = 0,
        P.uncorr.samp = 0,
        cpi = NA_real_,
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
  if(is.null(lengthValue)){
    P.boots = mapply(FUN = igr_prod.sample,
                     df = bootList,
                     # lengthValue = lengthValue,
                     massValue = massValue,
                     abunValue = abunValue,
                     dateCol = dateCol,
                     repCol = repCol,
                     wrap = wrap,
                     full = FALSE,
                     MoreArgs = list(dateDf = dateDf,
                                     growthDf = growthDf))

  } else{
  P.boots = mapply(FUN = igr_prod.sample,
                   df = bootList,
                   lengthValue = lengthValue,
                   massValue = massValue,
                   abunValue = abunValue,
                   dateCol = dateCol,
                   repCol = repCol,
                   wrap = wrap,
                   full = FALSE,
                   MoreArgs = list(dateDf = dateDf,
                                   growthDf = growthDf))
  }

  #### create SAMPLE information to export as summary ####
  sampSummary = create_sample_summary(df = taxaSampleListMass,
                                      wrap = wrap,
                                      abunValue = abunValue,
                                      massValue = massValue,
                                      dateCol = dateCol,
                                      repCol = repCol,
                                      ...)
  sampSummary <- merge(sampSummary, P.samp$dfDateAgg, by = dateCol, all.x = TRUE)
  #estimate the sample PB
  pb = P.samp$P.ann.samp/P.samp$B.ann.mean

  if(taxaSummary == FALSE){
    taxaSummary <- NULL
  } else if (taxaSummary == TRUE) {
    # # create a list for output
    taxaSummary <- list(
      taxonID = taxaInfo$taxonID,
      method = "igr",
      P.ann.samp = P.samp$P.ann.samp,
      pb = pb,
      meanN = P.samp$N.ann.mean,
      meanB = P.samp$B.ann.mean,
      meanIndMass = P.samp$B.ann.mean/ P.samp$N.ann.mean,
      datesInfo = sampSummary
    )
  }
  return(assign(speciesName, list(P.boots = P.boots,
              taxaSummary = taxaSummary)))

}
