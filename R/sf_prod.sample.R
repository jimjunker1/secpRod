#' @name sf_prod.sample
#' @title sf_prod.sample
#' @description This function calculates taxa production based on the size-frequency method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the size class including lengthClass, massClass, bin_min, bin_max, midpoint
#' @param massValue character string identifying the column name of the mass value
#' @param abunValue character string identifying the column name of the density value
#' @param dateCol character string identifying the column name of the sample date information. This is distinguished from *Value parameters in that Values may be used to maintain units provenance in the future. This may also change.
#' @param repCol character string identifying the column name of the replicate information.
#' @param wrap logical. should the calculations be wrapped by adding an additional date to make a full year?
#' @param cpi integer. The cohort production interval.
#' @param full logical. should the full summary be returned with mean and sd
#' @param ... additional arguments passed to function
#' @return list object with annual production, mean biomass, and mean abundance
#' @importFrom stats filter
#' @importFrom stats aggregate
#' @export
sf_prod.sample <- function(df = NULL,
                           sizesDf = NULL,
                           lengthValue = NULL,
                           massValue = 'afdm_mg',
                           abunValue = 'density',
                           dateCol = 'dateID',
                           repCol = 'repID',
                           wrap = FALSE,
                           cpi = NULL,
                           full = TRUE,
                           ...) {

  #### tests ####
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates
  df[["biomass"]] <- df[[abunValue]] * df[[massValue]]
  N.ann.list = estimate_ann_stats(df = df,
                                  var = abunValue,
                                  massValue = massValue,
                                  abunValue = abunValue,
                                  dateCol = dateCol,
                                  repCol = repCol,
                                  wrap = wrap)
  B.ann.list = estimate_ann_stats(df = df,
                                  var = "biomass",
                                  massValue = massValue,
                                  abunValue = abunValue,
                                  dateCol = dateCol,
                                  repCol = repCol,
                                  wrap = wrap)
 # df[["biomass_mean"]] <- NULL
 if(B.ann.list[["biomass_mean"]] == 0){
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

   df$lengthClass <- NA_real_
   df$massClass <- NA_real_
 if(is.null(lengthValue)){

   if(all(is.na(sizesDf$bin_min))){
     df$massClass <- df[[massValue]]
   }else{
        # x  <- sapply(df[[massValue]], FUN = function(a) sizesDf[apply(sizesDf,1,FUN = function(b) between(a,b[['bin_min']], b[['bin_max']])),'massClass'])
        breaks <- c(sizesDf$bin_min, tail(sizesDf$bin_max, 1))
        df$massClass <- sizesDf$massClass[findInterval(df[[massValue]], breaks, rightmost.closed = TRUE)]
        df$massClass[is.na(df$massClass)] <- NA
   }
 } else{
   df$lengthClass <- df[[lengthValue]]
   df$massClass <- df[[massValue]]
 }
   classedDf <- c()
   if(all(is.na(df$lengthClass))){
     densityRepAggForm <- paste0(abunValue,"~taxonID+massClass+",repCol)
   } else{
     densityRepAggForm <- paste0(abunValue,"~taxonID+lengthClass+massClass+",repCol)
   }
   repClassedDf <- stats::aggregate(formula(densityRepAggForm), data = df, FUN = sum, na.action = na.omit)
   densityAggForm <- gsub(paste0("\\+",repCol), "", densityRepAggForm)
   classedDf[[abunValue]] <- stats::aggregate(formula(densityAggForm), data = repClassedDf, FUN = mean)[[abunValue]]

   if(all(is.na(df$lengthClass))){
     massAggForm <- paste0(massValue,"~taxonID+massClass")
   } else{
     massAggForm <- paste0(massValue,"~taxonID+lengthClass+massClass")
   }
   # classedDf[[massValue]] <- stats::aggregate(formula())
  #### calculate SAMPLE annual production ####
  # Create a matrix with these 8 columns:
  # [1] size class (mm or mass),
  # [2] mean density for all samples throughout year (number m^-2),
  # [3] individual mass (mg AFDM) if mass is used for size class this will be duplicate of [1],
  # [4] density loss
  # [5] biomass density for each size class (rows),
  # [6] mass at loss (mean mass between size classes),
  # [7] biomass loss,
  # [8] multiply by # of size classes

  # 0) create the matrix object
  sfTab <- data.frame(matrix(NA_real_, nrow(sizesDf)+1, 9))
  names(sfTab) <- c('lengthClass','massClass', abunValue, 'ind.mass', paste0(abunValue,'.loss'),'biomass','mass.at.loss','biomass.loss','multiply.by.no.sizes')

  # 1) & 2) add lengthClass and massClass columns
  sfTab[,1] <- c(sizesDf[[1]],NA)
  sfTab[,2] <- c(sizesDf[[2]],NA)

  # 3) add in the density column
  sfTab[,3] <- c(classedDf[[abunValue]],NA)

  # 4) add in ind.mass column
  sfTab[,4] <- c(sizesDf[[2]],NA)

  # 5) add in the density.loss column
  sfTab[,5] <- c(NA, (diff(sfTab[,3]) *-1))
  sfTab[nrow(sfTab),5] <- sfTab[(nrow(sfTab)-1),3]

  # 6) add in the biomass column
  sfTab[,6] <- sfTab[['ind.mass']] * sfTab[[abunValue]]

  # 7) add in mass at loss column
  sfTab[,7] <- stats::filter(sfTab[['ind.mass']], c(1,1)/2, sides = 1)
  sfTab[nrow(sfTab),7] <- sfTab[(nrow(sfTab)-1),4]

  # 8) add in the biomass loss
  sfTab[,8] <- sfTab[['mass.at.loss']] * sfTab[[paste0(abunValue,'.loss')]]

  # 9) multiply by # of size classes
  sfTab[,9] <- sfTab[["biomass.loss"]] * (nrow(sizesDf)-1)
  # if the first biomass.loss x size classes column is negative, set to 0 as suggested by Benke, Huryn, and Junker 2026
  if(sfTab[2,9] < 0) sfTab[2,9] <- 0

  # Set negative values to zero only if no positive values precede them and they occur below a non-positive value (i.e., negative or zero) in the column of biomass * number of size classes (mg AFDM m^-2) as Benke & Huryn (2007) suggest
  # for (s in 2:dim(SF.int)[1]) {
  #   if (SF.int[s, 4] < 0 & sum(SF.int[1:(s - 1), 4] > 0) == 0) {
  #     SF.int[s, 4] <- 0
  #   }
  # }
  # Calculate "uncorrected" production by summing all values in the the column of biomass * number of size classes
  P.uncorr.samp <- sum(sfTab[, 9], na.rm = TRUE)
  # Calculate annual production using the cohort production interval (cpi) given in days for this taxon
  P.ann.samp <- P.uncorr.samp * (365 / cpi)

}
  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,
                P.uncorr.samp = P.uncorr.samp,
                B.ann.mean = B.ann.list[["biomass_mean"]],
                B.ann.sd = B.ann.list[["biomass_sd"]],
                N.ann.mean = N.ann.list[[paste0(abunValue,"_mean")]],
                N.ann.sd = N.ann.list[[paste0(abunValue,"_sd")]]))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                P.uncorr.samp = P.uncorr.samp,
                B.ann.samp = B.ann.list[["biomass_mean"]],
                N.ann.samp = N.ann.list[[paste0(abunValue,"_mean")]]))
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
