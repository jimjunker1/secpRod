#' @name rs_prod.sample
#' @title rs_prod.sample
#' @description This function calculates taxa production based on the removal-summation method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the
#' @param massValue character string identifying the column name of the mass value
#' @param massLabel character string identifying the column name of the mass value
#' @param full logical. should the full summary be returned with mean and sd
#' @param ... additional arguments passed to function
#' @return list object with taxa summary of the sampled data
#' @export
rs_prod.sample <- function(df = NULL,
                           sizesDf = NULL,
                           massValue = NULL,
                           massLabel = NULL,
                           full = TRUE,
                           ...) {

  #### tests ####
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates
  df[[massLabel]] <- df$n_m2 * df[[massValue]]
  N.ann.list = estimate_ann_stats(df,
                     var = 'n_m2')
  B.ann.list = estimate_ann_stats(df,
                                  var = massLabel)
 df[[massLabel]] <- NULL
 if(B.ann.list[[eval(paste0(massLabel,"_mean"))]] == 0){
   if(full == TRUE){
     return(list(P.ann.samp = 0,
                 P.uncorr.samp = 0,
                 B.ann.mean = 0,
                 B.ann.sd = NA_real_,
                 N.ann.mean = 0,
                 N.ann.sd = NA_real_))
   } else{
     return(list(P.ann.samp = 0,
                 B.ann.samp = 0,
                 N.ann.samp = 0))
   }

 } else{


  #### calculate SAMPLE annual production ####
  # Create a matrix with these 4 columns: individual length (mm), mean density for all samples throughout year (number m^-2), individual mass (mg AFDM), and biomass (mg AFDM m^-2) for each size class (rows)
  SF <- matrix(0, length(unique(unlist(df$lengthClass))), 4)

  SF[, c(1,3)] <- c(
    sizesDf[[1]], # lengthClass
    sizesDf[[2]] # massClass
  )
  SF[, 2] <- unname(unlist(aggregate(df$n_m2, by = list(df$lengthClass), mean, na.rm = TRUE)[2]))
  SF[, 4] <- SF[, 2] * SF[, 3]

  # Create a matrix with these 4 columns: number lost (number m^-2), individual mass at loss (mg AFDM), biomass lost (mg AFDM m^-2), and biomass lost * number size classes (mg AFDM m^-2) for each transistion between size classes (rows)
  SF.int <- matrix(0, length(unique(unlist(df$lengthClass))), 4)
  # Calculate the number lost between size classes, but subtract zero from the mean number in the largest size class for the "final" transition out of the largest size class
  SF.int[, 1] <- c(-diff(SF[, 2]), (SF[dim(SF)[1], 2] - 0))
  # Calculate the geometric mean of individual masses between size classes, but use the individual mass of the largest size class for the "final" transition out of the largest size class, as Benke & Huryn (2007) suggest
  SF.int[, 2] <- c((SF[(1:(dim(SF)[1] - 1)), 3] * SF[(2:dim(SF)[1]), 3])^(1 / 2), SF[dim(SF)[1], 3])
  SF.int[, 3] <- SF.int[, 1] * SF.int[, 2]
  SF.int[, 4] <- SF.int[, 3] * max(sizesDf[[1]])
  # If the first value in the column of biomass * number of size classes (mg AFDM m^-2) is negative, set it to zero
  if (SF.int[1, 4] < 0) {
    SF.int[1, 4] <- 0
  }
  # Set negative values to zero only if no positive values precede them and they occur below a non-positive value (i.e., negative or zero) in the column of biomass * number of size classes (mg AFDM m^-2) as Benke & Huryn (2007) suggest
  for (s in 2:dim(SF.int)[1]) {
    if (SF.int[s, 4] < 0 & sum(SF.int[1:(s - 1), 4] > 0) == 0) {
      SF.int[s, 4] <- 0
    }
  }

  # Calculate "uncorrected" production by summing all values in the the column of biomass * number of size classes (mg AFDM m^-2)
  P.uncorr.samp <- sum(SF.int[, 4])

}
  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,
                P.uncorr.samp = P.uncorr.samp,
                B.ann.mean = B.ann.list[[paste0(massLabel,"_mean")]],
                B.ann.sd = B.ann.list[[paste0(massLabel,"_sd")]],
                N.ann.mean = N.ann.list$n_m2_mean,
                N.ann.sd = N.ann.list$n_m2_sd))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                B.ann.samp = B.ann.list[[paste0(massLabel,"_mean")]],
                N.ann.samp = N.ann.list$n_m2_mean))
  }
  # #### create SAMPLE information to export as summary ####
  # # summarise sample sizes across dates
  # sampDatesInfo <- setNames(unique(aggregate(df[c("repID")], by = list(df$dateID, df$lengthClass), count)[c(1, 3)]), c("dateID", "N"))
  # if (wrap) {
  #   temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
  #   temp[["N"]] <- NA
  #   sampDatesInfo <- rbind(sampDatesInfo, temp)
  # }
  # # summarise the sample abundance N across all dates and size classes
  # Nmean <- setNames(cleanAggDf(aggregate(df, by = list(df$dateID, df$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", "n_m2_mean"))
  # Nsd <- setNames(cleanAggDf(aggregate(df["n_m2"], by = list(df$dateID, df$lengthClass), sd, na.rm = TRUE)), nm = "n_m2_sd")
  # Nbind <- cbind(Nmean, Nsd)
  # Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  # NmeanTab <- as.data.frame.matrix(xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  # NsdTab <- as.data.frame.matrix(xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  # NdatesInfo <- setNames(aggregate(Nmean["n_m2_mean"], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
  # # if wrap equals true create another
  # if (wrap) {
  #   temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
  #   temp[["n_m2_mean"]] <- mean(c(NdatesInfo[1, "n_m2_mean"], NdatesInfo[nrow(NdatesInfo), "n_m2_mean"]))
  #   NdatesInfo <- rbind(NdatesInfo, temp)
  # }
  # # summarise the sample biomasses across all dates and size classes
  # # create the sizeclass biomass for all rows
  # df[[massLabel]] <- unlist(df[, "n_m2"]) * unlist(df[, massValue])
  # # do the aggregating
  # Bmean <- setNames(cleanAggDf(aggregate(df[c("dateID", "lengthClass", massLabel)], by = list(df$dateID, df$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  # Bsd <- setNames(cleanAggDf(aggregate(df[massLabel], by = list(df$dateID, df$lengthClass), sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  # Bbind <- cbind(Bmean, Bsd)
  # Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  # meanBform <- as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  # sdBform <- as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  # BmeanTab <- as.data.frame.matrix(xtabs(meanBform, Bbind))
  # BsdTab <- as.data.frame.matrix(xtabs(sdBform, Bbind))
  # BdatesInfo <- setNames(aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # # if wrap equals true
  # if (wrap) {
  #   temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
  #   temp[[eval(paste0(massLabel, "_mean"))]] <- mean(c(BdatesInfo[1, eval(paste0(massLabel, "_mean"))], BdatesInfo[nrow(BdatesInfo), eval(paste0(massLabel, "_mean"))]))
  #
  #   BdatesInfo <- rbind(BdatesInfo, temp)
  # }
  # # create the full summary
  # datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, NdatesInfo, BdatesInfo))
  #
  # if(taxaSummary == "none"){
  #
  # } else if (taxaSummary == "full") {
  #   # # create a list for output
  #   taxaSummary <- list(
  #     summaryType = "full",
  #     taxonID = taxaInfo$taxonID,
  #     method = "sf",
  #     P.ann.samp = P.ann.samp$P.ann.samp,
  #     P.uncorr.samp = P.uncorr.samp$P.uncorr.samp,
  #     cpi = cpi,
  #     meanN = mean(unlist(datesInfo$n_m2_mean)),
  #     meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
  #     meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
  #     Nmean = NmeanTab,
  #     Nsd = NsdTab,
  #     Bmean = BmeanTab,
  #     Bsd = BsdTab,
  #     datesInfo = datesInfo
  #   )
  # } else if(taxaSummary == "short"){
  #   taxaSummary <- list(
  #     summaryType = "short",
  #     taxonID = taxaInfo$taxonID,
  #     method = "sf",
  #     P.ann.samp = P.ann.samp$P.ann.samp,
  #     cpi = cpi,
  #     meanN = mean(unlist(datesInfo$n_m2_mean)),
  #     meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
  #     meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
  #     datesInfo = datesInfo
  #   )
  # }
  # #   assign(taxaInfo$taxonID, list())
  # # # add taxainformation
  # #   assign(taxaInfo$taxonID,
  # #          within(eval(as.symbol(taxaInfo$taxonID)),{
  # #            taxaInfo <- taxaInfo
  # #            }
  # #            )
  # #          )
  #
  # return(taxaSummary)
}
