#' @name is_prod.sample
#' @title is_prod.sample
#' @description This function calculates taxa production based on the increment-summation method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the
#' @param massValue character string identifying the column name of the mass value
#' @param massLabel character string identifying the column name of the mass value
#' @param full logical. should the full summary be returned with mean and sd
#' @param ... additional arguments passed to function
#' @return list object with taxa summary of the sampled data
#' @importFrom stats filter
#' @importFrom stats aggregate
#' @importFrom zoo na.approx
#' @export
is_prod.sample <- function(df = NULL,
                           dateDf = dateDf,
                           massValue = 'afdm_mg',
                           abunValue = 'density',
                           dateCol = 'dateID',
                           repCol = 'repID',
                           wrap = FALSE,
                           full = TRUE,
                           ...) {
  #### tests ####

  #### end tests ####
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates for summaries
  df[["biomass"]] <- df[[abunValue]] * df[[massValue]]
  N.ann.list = estimate_ann_stats(df, var = abunValue,
                                  massValue = 'afdm_mg',
                                  abunValue = 'density',
                                  dateCol = 'dateID',
                                  repCol = 'repID',
                                  wrap = wrap)

  B.ann.list = estimate_ann_stats(df, var = "biomass",
                                  massValue = 'afdm_mg',
                                  abunValue = 'density',
                                  dateCol = 'dateID',
                                  repCol = 'repID',
                                  wrap = wrap)

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

  #### calculate SAMPLE annual production ####
  # Create a matrix with these 8 columns:
  # [1] date
  # [2] density,
  # [3] mean individual mass,
  # [4] biomass,
  # [5] gross individual growth,
  # [6] mean density,
  # [7] interval P,
  # [8] daily P,
  # [9] growth rate (d^-1),

  # 0) create the matrix object
  isTab <- data.frame(matrix(0, length(unique(unlist(df[[dateCol]]))), 9))
  names(isTab) <- c(dateCol, abunValue, 'ind.mass','biomass','mean.growth','density.mean','p.int','p.daily','g.daily')

  # 1) add date column
  isTab[,1] <- sort(unique(unlist(df[[dateCol]])))

  # 2) add density column
  ## first sum across all size classes within each replicate
  densityAgg = aggregate(formula(paste0(abunValue,"~",dateCol,"+",repCol)), data = df, FUN = sum, na.action = na.omit)
  ## then take the mean for each date
  isTab[,2] <- aggregate(formula(paste0(abunValue,"~",dateCol)), data = densityAgg, FUN = mean, na.action = na.omit, simplify = TRUE)[,2]

  # 3) add mean individual mass column
  ## first create a weighted mean across all size classes within a replicate
  ### create size class sum density
  dateMassSums <- aggregate(formula(paste0(abunValue,"~",massValue,"+",dateCol)), data = df, FUN = sum, na.action = na.omit)
  ### create date-level total density
  dateSums <- aggregate(formula(paste0(abunValue,"~",dateCol)), data = df, FUN = sum, na.action = na.omit)
  ### merge size class sum density and date-level total density
  dateMerge <- merge(dateMassSums, dateSums, by = eval(dateCol))
  ### calculate the weights for each size class based on relative density
  dateMerge[['weights']] <- dateMerge$density.x/dateMerge$density.y
  ### weight the size class by relative weights
  dateMerge[['w.mass']] <- dateMerge$afdm_mg * dateMerge$weights
  ### sum the weighted size class within sampling date
  massAgg <- setNames(aggregate(formula(paste0('w.mass ~',dateCol)), data = dateMerge, FUN = sum, na.action = na.omit), nm = c(dateCol, 'w.mass'))
  massAgg <- merge(isTab[dateCol], massAgg, by = dateCol, all.x = TRUE)
  naMasses <- which(is.na(massAgg$w.mass))
  if(length(naMasses) > 0){
    if(all(length(naMasses) == 1 & naMasses == length(massAgg$w.mass))){
      massAgg$w.mass[length(massAgg$w.mass)] <- massAgg$w.mass[length(massAgg$w.mass)-1]
    } else if(any(diff(naMasses) == 1)){
      if(any(length(massAgg$w.mass) %in% naMasses)){
        massAgg$w.mass[length(massAgg$w.mass)] <- massAgg$w.mass[length(massAgg$w.mass)-1] <- massAgg$w.mass[length(massAgg$w.mass)-2]
      } else if(length(massAgg$w.mass) %ni% naMasses){
        naLocs = which(diff(naMasses) == 1)
        massAgg$w.mass[c((naLocs-1), naLocs)] <- zoo::na.approx(c(naLocs-2, naLocs -1, naLocs, naLocs +1))
      }} else{
  for(i in 1:length(naMasses)){
    massAgg$w.mass[naMasses[i]] <- (massAgg$w.mass[naMasses[i]-1]+massAgg$w.mass[naMasses[i]+1])/2
  }
    }
    }
  ## set the weighted mean individual mass for each date
  isTab[,3] <- massAgg$w.mass

  # 4) calculate biomass for each date
  isTab[,4] <- isTab[,2] * isTab[,3]

  # 5) add the gross individual growth column
  ## create a lagged mean using the massAgg object
  isTab[,5] <- c(NA,diff(massAgg$w.mass))

  # 6) create a lagged mean using the density calculations
  isTab[,6] <- stats::filter(isTab[,2], c(1,1)/2, sides = 1)

  # 7) estimate interval production
  isTab[,7] <- isTab[,5] * isTab[,4]

  # 8) estimate the daily production rates by merging interval lengths
  isTab[,8] <- isTab[,7] / c(NA, dateDf$int_days)

  # 9) estimate growth rate by merging interval lengths
  isTab[,9] <- isTab[,5] / c(NA, dateDf$int_days)

  # Calculate secondary production by summing initial biomass and all production values in isTab[,7]
  P.ann.samp <- isTab[1,4] + sum(isTab[, 7], na.rm = TRUE)

}
  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,
                B.ann.mean = B.ann.list[["biomass_mean"]],
                B.ann.sd = B.ann.list[["biomass_sd"]],
                N.ann.mean = N.ann.list[[paste0(abunValue,"_mean")]],
                N.ann.sd = N.ann.list[[paste0(abunValue,"_sd")]]))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                B.ann.samp = B.ann.list[["biomass_mean"]],
                N.ann.samp = N.ann.list[[paste0(abunValue,"_mean")]]))
  }
}
