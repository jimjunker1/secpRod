#' @description
#' This function calculates secondary production with the size-frequency method.
#' @title calc_prod_sf
#' @param taxaSampleListMass description
#' @param taxaInfo data frame of taxonomic information for calculating production
#' @param bootNum integer. How many bootstrap samples should be constructed
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
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats formula
#' @export

calc_prod_sf <- function(taxaSampleListMass= NULL,
                         taxaInfo = NULL,
                         bootNum = NULL,
                         # dateDf = NULL,
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

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  taxaInfo = taxaInfo[which(taxaInfo$taxonID == speciesName),]
  ## function prep ##
  ### make a list of key variables to pass to sample function
  #### create the sizesDf object
  #### Create size class breakdown if lengthValue not provided, i.e. just masses are present ####
    if(is.null(lengthValue)){
      # make a guess if the masses are classes or continuous
      # determine the number of unique size classes relative to total density?
      numSizes = length(unique(taxaSampleListMass[[massValue]]))
      totalDensity = sum(taxaSampleListMass[[abunValue]], na.rm = TRUE)
      if(any(numSizes/totalDensity > 0.2 | numSizes > 30)){
      # If continuous, bin the masses into X bins
      # X should be as large as possible, but not too big that there are bins with 0s
      sizes = as.numeric(taxaSampleListMass[[massValue]])
      # get min and max sizes
      min_size = min(sizes, na.rm = TRUE);max_size = max(sizes, na.rm = TRUE)
      # set breaks for equally spaced bins within this range with a minimum of 6

      zeroBins = sapply(15:6, FUN = function(x){
        breaks = seq(min_size, max_size, length.out = x +1)
        bins = cut(sizes, breaks = breaks, include.lowest = TRUE)
        zeros = !any(unlist(table(bins)) == 0)
        negative = all(diff(unname(table(bins))) < 0)
        return(data.frame(zeros = list(zeros),
                          negative = list(negative)))
        })
      suppressWarnings(
        # if all bins have either zeros or non-negative changes
        # set at minimum of 6. Warnings are suppressed for coercing list to logical
        if(all((15:6)[apply(zeroBins,2, function(x) all(x))])){
        b = 6
      } else{
        b = (15:6)[which.max(apply(zeroBins,2, function(x) all(x)))]
      }
      )

      breaks = seq(min_size, max_size, length.out = b+1)
      bins = cut(sizes, breaks = breaks, include.lowest = TRUE)
      sizesDf = get_bin_widths(levels(bins))
      sizesDf$lengthClass <- NA
      sizesDf$massClass <- as.numeric(unlist(sizesDf$midpoint))
      sizesDf = sizesDf[,c('lengthClass','massClass','bin_min','bin_max')]

      } else{ # if the data appear to be already in classes
      sizesDf$lengthClass <- sizesDf$bin_min <- sizesDf$bin_max <- NA
      sizesDf$massClass = unique(taxaSampleListMass[, c(eval(massValue))])
      sizesDf = sizesDf[,c('lengthClass','massClass','bin_min','bin_max')]

      }
    } else{

      sizesDf = unique(taxaSampleListMass[, c(eval(lengthValue), eval(massValue))])
      sizesDf$bin_min <- sizesDf$bin_max <- NA
      sizesDf = sizesDf[,c('lengthClass','massClass','bin_min','bin_max')]

    }

  ## set the funcList to pass to sample production function
  funcList = list(
    df = taxaSampleListMass,
    sizesDf = sizesDf,
    massValue = massValue,
    abunValue = abunValue,
    wrap = wrap
  )
  # calculate the production from the observed samples
  taxaCPI <- mean(c(taxaInfo$min.cpi, taxaInfo$max.cpi))
  funcList = c(funcList, list(cpi = taxaCPI))
  P.samp = do.call(sf_prod.sample, args = funcList)
  if(P.samp$P.ann.samp == 0){
    if(taxaSummary == FALSE){
      taxaSummary <- NULL
    } else if (taxaSummary == TRUE) {
      # # create a list for output
      taxaSummary <- list(
        taxonID = taxaInfo$taxonID,
        method = "sf",
        P.ann.samp = 0,
        P.uncorr.samp = 0,
        cpi = NA_real_,
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

  cpiBoots = sample(seq.int(from = as.integer(taxaInfo$min.cpi), to = as.integer(taxaInfo$max.cpi), by = 1), as.integer(bootNum), replace = TRUE)


  # P.boots = vector('list', length = bootNum)
  # for(i in 1:length(bootList)){
  #   P.boots[[i]] = sf_prod.sample(
  #                  df = bootList[[i]],
  #                  sizesDf = sizesDf,
  #                  lengthValue = lengthValue,
  #                  massValue = massValue,
  #                  abunValue = abunValue,
  #                  dateCol = dateCol,
  #                  repCol = repCol,
  #                  cpi = cpiBoots[i],
  #                  wrap = wrap,
  #                  full = FALSE,
  #                  )
  # }
  # debug(sf_prod.sample)
  P.boots = mapply(FUN = sf_prod.sample,
                   df = bootList,
                   # lengthValue = lengthValue,
                   massValue = massValue,
                   abunValue = abunValue,
                   dateCol = dateCol,
                   repCol = repCol,
                   cpi = cpiBoots,
                   wrap = wrap,
                   full = FALSE,
                   MoreArgs = list(sizesDf = sizesDf))

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
      method = "sf",
      P.ann.samp = P.samp$P.ann.samp,
      P.uncorr.samp = P.samp$P.uncorr.samp,
      cpi = taxaCPI,
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
