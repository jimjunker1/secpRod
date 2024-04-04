#' @name sf_prod.boot
#' @description This function calculates taxa production based on the size-frequency method
#' @param taxaSampleListMass a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param taxaInfo a data.frame of the important taxonomic information for converting length to mass, growth rate formulas, cohort production intervals (cpi), etc.
#' @param fullTaxaSummary logical. Should the taxa summary returned be full or abridged
#' @return list object with taxa summary of the sampled data
#' @export sf_prod.boot
sf_prod.boot <- function(taxaSampleListMass = NULL, taxaInfo = NULL,...) {

  #### tests ####
  if (is.null(taxaInfo)) stop("`taxaInfo` must be provided. Currently is NULL")

  #### function prep ####
  ## Grab a list of key variables
  massValue <- rev(names(taxaSampleListMass))[1]
  massLabel <- paste0(massValue, "_m2")
  # wrap dates
  if(wrap)
  dateDf <- wrap_dates(taxaSampleListMass, wrap = wrap,...)

  # dataframe of sizes and masses

  # fullSizesDf = convert_length_to_mass(taxaSampleListMass, infoCols, taxaInfo)
  sizesDf <- unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])])

  # allocate vectors
  # cpi vector
  cpiVec = vector(mode = 'numeric', length = bootNum)

cpiVec = runif(bootNum, min = taxaInfo$min.cpi, max = taxaInfo$max.cpi)

# # Remove any stochasticity when calculating SAMPLE production by using the mean of the min & max cpi values
# cpi <- mean(c(taxaInfo$min.cpi, taxaInfo$max.cpi))
#
# calc_prod_sf(df = taxaSampleListMass, cpi = cpi)

  #### create SAMPLE information to export as summary ####
  # summarise sample sizes across dates
  sampDatesInfo <- setNames(unique(aggregate(taxaSampleListMass[c("repID")], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), count)[c(1, 3)]), c("dateID", "N"))
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[["N"]] <- NA
    sampDatesInfo <- rbind(sampDatesInfo, temp)
  }
  # summarise the sample abundance N across all dates and size classes
  Nmean <- setNames(cleanAggDf(aggregate(taxaSampleListMass, by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", "n_m2_mean"))
  Nsd <- setNames(cleanAggDf(aggregate(taxaSampleListMass["n_m2"], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), sd, na.rm = TRUE)), nm = "n_m2_sd")
  Nbind <- cbind(Nmean, Nsd)
  Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  NmeanTab <- as.data.frame.matrix(xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  NsdTab <- as.data.frame.matrix(xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  NdatesInfo <- setNames(aggregate(Nmean["n_m2_mean"], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
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
  Bmean <- setNames(cleanAggDf(aggregate(taxaSampleListMass[c("dateID", "lengthClass", massLabel)], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  Bsd <- setNames(cleanAggDf(aggregate(taxaSampleListMass[massLabel], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  Bbind <- cbind(Bmean, Bsd)
  Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  meanBform <- as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  sdBform <- as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  BmeanTab <- as.data.frame.matrix(xtabs(meanBform, Bbind))
  BsdTab <- as.data.frame.matrix(xtabs(sdBform, Bbind))
  BdatesInfo <- setNames(aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # if wrap equals true
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[[eval(paste0(massLabel, "_mean"))]] <- mean(c(BdatesInfo[1, eval(paste0(massLabel, "_mean"))], BdatesInfo[nrow(BdatesInfo), eval(paste0(massLabel, "_mean"))]))

    BdatesInfo <- rbind(BdatesInfo, temp)
  }
  # create the full summary
  datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, NdatesInfo, BdatesInfo))

  # # create a list for output
  if (fullTaxaSummary) {
    taxaSummary <- list(
      summaryType = "full",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.ann.samp,
      P.uncorr.samp = P.uncorr.samp,
      cpi = cpi,
      meanN = mean(unlist(datesInfo$n_m2_mean)),
      meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
      meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
      Nmean = NmeanTab,
      Nsd = NsdTab,
      Bmean = BmeanTab,
      Bsd = BsdTab,
      datesInfo = datesInfo
    )
  } else {
    taxaSummary <- list(
      summaryType = "short",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.ann.samp,
      cpi = cpi,
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

  return(taxaSummary)
}

