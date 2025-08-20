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
                         bootList = NULL,...) {

  ## tests ##

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  # ## function prep ##
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    # sizesDf = unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])])
    # sizesDf = unique(taxaSampleListMass[, c(eval(lengthValue), eval(massValue))]),
    massValue = massValue,
    abunValue = abunValue,
    dateDf = dateDf,
    dateCol = dateCol# massLabel = massLabel
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
  P.boots = mapply(FUN = is_prod.sample,
                   df = bootList,
                   # sizesDf = lapply(1:bootNum, function(x) funcList$sizesDf),
                   massValue = massValue,
                   abunValue = abunValue,
                   # dateDf = dateDf,
                   # dateDf = lapply(1:bootNum, function(x) dateDf),
                   dateCol = dateCol,
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
  # # summarise sample sizes across dates
  # ## count the number of replicates across dates
  # sampDatesInfo <- stats::setNames(stats::aggregate(taxaSampleListMass[[repCol]], by = list(taxaSampleListMass[[dateCol]]), FUN = function(x) length(unique(x))), c("dateID", "N"))
  # if (wrap) {
  #   temp <- wrap_dates(df = taxaSampleListMass, dateCol = dateCol, wrapDate = TRUE)[dateCol]
  #   temp[["N"]] <- NA
  #   sampDatesInfo <- rbind(sampDatesInfo, temp[nrow(temp),])
  #   temp[["N"]] <- NULL
  # }
  # ## get sample abundance across dates dates
  # densityAgg <- aggregate(formula(paste0(abunValue,"~",dateCol,"+",repCol)), data = taxaSampleListMass, FUN = sum, na.rm = TRUE)
  # ## then take the mean and sd for each date
  # Nmean <- stats::setNames(aggregate(formula(paste0(abunValue,"~",dateCol)), data = densityAgg, FUN = mean, na.rm = TRUE), nm = c("dateID", paste0(abunValue,"_mean")))
  # Nsd <- stats::setNames(aggregate(formula(paste0(abunValue,"~",dateCol)), data = densityAgg, FUN = stats::sd, na.rm = TRUE), nm = c("dateID", paste0(abunValue,"_sd")))
  # # bind them together
  # Nbind <- merge(Nmean, Nsd, by = dateCol)
  # # Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  # # NmeanTab <- as.data.frame.matrix(stats::xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  # # NsdTab <- as.data.frame.matrix(stats::xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  # # NdatesInfo <- stats::setNames(stats::aggregate(Nmean[[paste0(abunValue,"_mean")]], by = list(Nmean[[dateCol]]), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
  # # if wrap equals true create another
  # if (wrap) {
  #   temp[[paste0(abunValue,"_mean")]] <- temp[[paste0(abunValue,"_sd")]] <- NA
  #   temp[nrow(temp), paste0(abunValue,"_mean")] <- (Nbind[1,paste0(abunValue,"_mean")] + Nbind[(nrow(Nbind)-1),paste0(abunValue,"_mean")])/2
  #   temp[nrow(temp), paste0(abunValue,"_sd")] <- sqrt((Nbind[1,paste0(abunValue,"_sd")]^2) + (Nbind[(nrow(Nbind)-1), paste0(abunValue,"_sd")]^2))
  #   Nbind <- rbind(Nbind, temp[nrow(temp),])
  #   temp[[paste0(abunValue,"_mean")]] <- temp[[paste0(abunValue,"_sd")]] <- NULL
  # }
  # # summarise the sample biomasses across all dates and size classes
  # # create the sizeclass biomass for all rows
  # taxaSampleListMass[["biomass"]] <- unlist(taxaSampleListMass[, abunValue]) * unlist(taxaSampleListMass[, massValue])
  # ## convert biomass NAs to 0
  # taxaSampleListMass[,"biomass"][is.na(taxaSampleListMass[,"biomass"])] <- 0
  # # do the aggregating
  # ## get sample abundance across dates dates
  # biomassAgg <- aggregate(formula(paste0("biomass~",dateCol,"+",repCol)), data = taxaSampleListMass, FUN = sum, na.rm = TRUE)
  # ## then take the mean and sd for each date
  # Bmean <- stats::setNames(aggregate(formula(paste0("biomass~",dateCol)), data = biomassAgg, FUN = mean, na.rm = TRUE), nm = c("dateID", "biomass_mean"))
  # Bsd <- stats::setNames(aggregate(formula(paste0("biomass~",dateCol)), data = biomassAgg, FUN = stats::sd, na.rm = TRUE), nm = c("dateID", "biomass_sd"))
  # # bind them together
  # Bbind <- merge(Bmean, Bsd, by = dateCol)
  # # Bmean <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[c("dateID", "lengthClass", massLabel)], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  # # Bsd <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[massLabel], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), stats::sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  # # Bbind <- cbind(Bmean, Bsd)
  # # Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  # # meanBform <- stats::as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  # # sdBform <- stats::as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  # # BmeanTab <- as.data.frame.matrix(stats::xtabs(meanBform, Bbind))
  # # BsdTab <- as.data.frame.matrix(stats::xtabs(sdBform, Bbind))
  # # BdatesInfo <- stats::setNames(stats::aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # # if wrap equals true
  # if (wrap) {
  #   temp[["biomass_mean"]] <- temp[["biomass_sd"]] <- NA
  #   temp[nrow(temp), "biomass_mean"] <- (Bbind[1,"biomass_mean"] + Bbind[(nrow(Bbind)-1),"biomass_mean"])/2
  #   temp[nrow(temp), "biomass_sd"] <- sqrt((Bbind[1,"biomass_sd"]^2) + (Bbind[(nrow(Bbind)-1), "biomass_sd"]^2))
  #   Bbind <- rbind(Bbind, temp[nrow(temp),])
  #   temp[["biomass_mean"]] <- temp[["biomass_sd"]] <- NULL
  # }
  # # create the full summary
  # datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, Nbind, Bbind))
  # #estimate the sample PB
  # pb = P.samp$P.ann.samp/mean(unlist(datesInfo[["biomass_mean"]]))
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
