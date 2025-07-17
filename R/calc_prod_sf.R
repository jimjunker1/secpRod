#' @description
#' This function calculates secondary production with the size-frequency method.
#' @title calc_prod_sf
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
#' @importFrom stats xtabs
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats as.formula
#' @export

calc_prod_sf <- function(taxaSampleListMass= NULL,
                         taxaInfo = NULL,
                         bootNum = NULL,
                         dateDf = NULL,
                         taxaSummary = 'full',
                         wrap = TRUE,
                         massValue = NULL,
                         massLabel = NULL,
                         bootList = NULL,...) {

  ## tests ##

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  taxaInfo = taxaInfo[which(taxaInfo$taxonID == speciesName),]
  # ## function prep ##
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    # sizesDf = unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])])
    sizesDf = unique(taxaSampleListMass[, c("lengthClass", eval(massValue))]),
    massValue = massValue,
    massLabel = massLabel
  )

  # calculate the production from the full samples
  taxaCPI <- mean(c(taxaInfo$min.cpi, taxaInfo$max.cpi))
  funcList = c(funcList, list(cpi = taxaCPI))
  P.samp = do.call(sf_prod.sample, args = funcList)
  if(P.samp$P.ann.samp == 0){
    if(taxaSummary == "none"){
      taxaSummary <- NULL
    } else if (taxaSummary == "full") {
      # # create a list for output
      taxaSummary <- list(
        summaryType = "full",
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
    } else if(taxaSummary == "short"){
      taxaSummary <- list(
        summaryType = "short",
        taxonID = taxaInfo$taxonID,
        method = "sf",
        P.ann.samp = 0,
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

  cpiBoots = sample(seq.int(from = as.integer(taxaInfo$min.cpi), to = as.integer(taxaInfo$max.cpi), by = 1), as.integer(bootNum), replace = TRUE)

  P.boots = mapply(FUN = sf_prod.sample,
                   df = bootList,
                   sizesDf = lapply(1:bootNum, function(x) funcList$sizesDf),
                   massValue = massValue,
                   massLabel = massLabel,
                   cpi = cpiBoots,
                   full = FALSE)
  #### create SAMPLE information to export as summary ####
  # summarise sample sizes across dates
  sampDatesInfo <- stats::setNames(unique(stats::aggregate(taxaSampleListMass[c("repID")], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), vers_count)[c(1, 3)]), c("dateID", "N"))
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[["N"]] <- NA
    sampDatesInfo <- rbind(sampDatesInfo, temp)
  }
  # summarise the sample abundance N across all dates and size classes
  Nmean <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass, by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", "n_m2_mean"))
  Nsd <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass["n_m2"], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), stats::sd, na.rm = TRUE)), nm = "n_m2_sd")
  Nbind <- cbind(Nmean, Nsd)
  Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  NmeanTab <- as.data.frame.matrix(stats::xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  NsdTab <- as.data.frame.matrix(stats::xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  NdatesInfo <- stats::setNames(stats::aggregate(Nmean["n_m2_mean"], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
  # if wrap equals true create another
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[["n_m2_mean"]] <- mean(c(NdatesInfo[1, "n_m2_mean"], NdatesInfo[nrow(NdatesInfo), "n_m2_mean"]))
    NdatesInfo <- rbind(NdatesInfo, temp)
  }
  # summarise the sample biomasses across all dates and size classes
  # create the sizeclass biomass for all rows
  taxaSampleListMass[[massLabel]] <- unlist(taxaSampleListMass[, "n_m2"]) * unlist(taxaSampleListMass[, massValue])
  # do the aggregating
  Bmean <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[c("dateID", "lengthClass", massLabel)], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  Bsd <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[massLabel], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), stats::sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  Bbind <- cbind(Bmean, Bsd)
  Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  meanBform <- stats::as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  sdBform <- stats::as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  BmeanTab <- as.data.frame.matrix(stats::xtabs(meanBform, Bbind))
  BsdTab <- as.data.frame.matrix(stats::xtabs(sdBform, Bbind))
  BdatesInfo <- stats::setNames(stats::aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # if wrap equals true
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[[eval(paste0(massLabel, "_mean"))]] <- mean(c(BdatesInfo[1, eval(paste0(massLabel, "_mean"))], BdatesInfo[nrow(BdatesInfo), eval(paste0(massLabel, "_mean"))]))

    BdatesInfo <- rbind(BdatesInfo, temp)
  }
  # create the full summary
  datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, NdatesInfo, BdatesInfo))
  #estimate the sample PB
  pb = P.samp$P.ann.samp/mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]]))
  if(taxaSummary == "none"){

  } else if (taxaSummary == "full") {
    # # create a list for output
    taxaSummary <- list(
      summaryType = "full",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.samp$P.ann.samp,
      P.uncorr.samp = P.samp$P.uncorr.samp,
      cpi = taxaCPI,
      pb = pb,
      meanN = mean(unlist(datesInfo$n_m2_mean)),
      meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
      meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
      Nmean = NmeanTab,
      Nsd = NsdTab,
      Bmean = BmeanTab,
      Bsd = BsdTab,
      datesInfo = datesInfo
    )
  } else if(taxaSummary == "short"){
    taxaSummary <- list(
      summaryType = "short",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.samp$P.ann.samp,
      cpi = taxaCPI,
      pb = pb,
      meanN = mean(unlist(datesInfo$n_m2_mean)),
      meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
      meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
      datesInfo = datesInfo
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
