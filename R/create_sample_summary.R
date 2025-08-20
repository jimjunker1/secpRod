#' @title create_sample_summary
#' @description A function to convert species-specific lengths to mass based on a user-provided length-mass equation form and variable values
#' @param df data.frame in long format for a single taxa. The data.frame should contain a species identifier column `taxonID` and a column of length bin categories `lengthClass`. `lengthClass` values must be numeric or coercible.
#' @param wrap logical (or coercible) indicating whether the production estimate should wrap the first and last dates to create a full annual cycle. If TRUE, this will create an additional sampling interval using the mean densities and masses to create a full annual data set.
#' @param abunValue description
#' @param massValue string of the column name containing the mass measurement
#' @param abunValue string of the column name containing the abundance or density measurement
#' @param dateCol string of the column name containing the date information. This can be either a recognized date object (e.g., Date, POSIX)
#' @param repCol string of the column name containing the replicate information
#' @param ... additional arguments passed to function
#' @returns a summary data.frame of the sample summary for N: sample size, density: mean and sd, biomass: mean and sd
#' @importFrom stats aggregate
#' @importFrom stats setNames
#' @importFrom stats formula
#' @export

create_sample_summary <- function(df = NULL, wrap = FALSE, abunValue = NULL, massValue = NULL, dateCol = NULL, repCol = NULL,...){
  #### create SAMPLE information to export as summary ####
  # summarise sample sizes across dates
  ## count the number of replicates across dates
  sampDatesInfo <- stats::setNames(stats::aggregate(df[[repCol]], by = list(df[[dateCol]]), FUN = function(x) length(unique(x))), c("dateID", "N"))
  if (wrap) {
    temp <- wrap_dates(df = df, dateCol = dateCol, wrapDate = TRUE)[dateCol]
    temp[["N"]] <- NA
    sampDatesInfo <- rbind(sampDatesInfo, temp[nrow(temp),])
    temp[["N"]] <- NULL
  }
  ## get sample abundance across dates dates
  densityAgg <- stats::aggregate(formula(paste0(abunValue,"~",dateCol,"+",repCol)), data = df, FUN = sum, na.rm = TRUE)
  ## then take the mean and sd for each date
  Nmean <- stats::setNames(stats::aggregate(formula(paste0(abunValue,"~",dateCol)), data = densityAgg, FUN = mean, na.rm = TRUE), nm = c("dateID", paste0(abunValue,"_mean")))
  Nsd <- stats::setNames(stats::aggregate(formula(paste0(abunValue,"~",dateCol)), data = densityAgg, FUN = stats::sd, na.rm = TRUE), nm = c("dateID", paste0(abunValue,"_sd")))
  # bind them together
  Nbind <- merge(Nmean, Nsd, by = dateCol)
  # Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  # NmeanTab <- as.data.frame.matrix(stats::xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  # NsdTab <- as.data.frame.matrix(stats::xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  # NdatesInfo <- stats::setNames(stats::aggregate(Nmean[[paste0(abunValue,"_mean")]], by = list(Nmean[[dateCol]]), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
  # if wrap equals true create another
  if (wrap) {
    temp[[paste0(abunValue,"_mean")]] <- temp[[paste0(abunValue,"_sd")]] <- NA
    temp[nrow(temp), paste0(abunValue,"_mean")] <- (Nbind[1,paste0(abunValue,"_mean")] + Nbind[(nrow(Nbind)-1),paste0(abunValue,"_mean")])/2
    temp[nrow(temp), paste0(abunValue,"_sd")] <- sqrt((Nbind[1,paste0(abunValue,"_sd")]^2) + (Nbind[(nrow(Nbind)-1), paste0(abunValue,"_sd")]^2))
    Nbind <- rbind(Nbind, temp[nrow(temp),])
    temp[[paste0(abunValue,"_mean")]] <- temp[[paste0(abunValue,"_sd")]] <- NULL
  }
  # summarise the sample biomasses across all dates and size classes
  # create the sizeclass biomass for all rows
  df[["biomass"]] <- unlist(df[, abunValue]) * unlist(df[, massValue])
  ## convert biomass NAs to 0
  df[,"biomass"][is.na(df[,"biomass"])] <- 0
  # do the aggregating
  ## get sample abundance across dates dates
  biomassAgg <- aggregate(formula(paste0("biomass~",dateCol,"+",repCol)), data = df, FUN = sum, na.rm = TRUE)
  ## then take the mean and sd for each date
  Bmean <- stats::setNames(aggregate(formula(paste0("biomass~",dateCol)), data = biomassAgg, FUN = mean, na.rm = TRUE), nm = c("dateID", "biomass_mean"))
  Bsd <- stats::setNames(aggregate(formula(paste0("biomass~",dateCol)), data = biomassAgg, FUN = stats::sd, na.rm = TRUE), nm = c("dateID", "biomass_sd"))
  # bind them together
  Bbind <- merge(Bmean, Bsd, by = dateCol)
  # Bmean <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[c("dateID", "lengthClass", massLabel)], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  # Bsd <- stats::setNames(cleanAggDf(stats::aggregate(taxaSampleListMass[massLabel], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), stats::sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  # Bbind <- cbind(Bmean, Bsd)
  # Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  # meanBform <- stats::as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  # sdBform <- stats::as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  # BmeanTab <- as.data.frame.matrix(stats::xtabs(meanBform, Bbind))
  # BsdTab <- as.data.frame.matrix(stats::xtabs(sdBform, Bbind))
  # BdatesInfo <- stats::setNames(stats::aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # if wrap equals true
  if (wrap) {
    temp[["biomass_mean"]] <- temp[["biomass_sd"]] <- NA
    temp[nrow(temp), "biomass_mean"] <- (Bbind[1,"biomass_mean"] + Bbind[(nrow(Bbind)-1),"biomass_mean"])/2
    temp[nrow(temp), "biomass_sd"] <- sqrt((Bbind[1,"biomass_sd"]^2) + (Bbind[(nrow(Bbind)-1), "biomass_sd"]^2))
    Bbind <- rbind(Bbind, temp[nrow(temp),])
    temp[["biomass_mean"]] <- temp[["biomass_sd"]] <- NULL
  }
  # create the full summary
  datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, Nbind, Bbind))

  return(datesInfo)

}
