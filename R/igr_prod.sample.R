#' @title igr_prod.sample
#' @description This function calculates taxa production based on the instantaneous growth method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the
#' @param full logical. should the full summary be returned with mean and sd
#' @param cpi integer. The cohort production interval.
#' @param ... additional arguments passed to function
#' @return list object with taxa summary of the sampled data
#' @export
igr_prod.sample <- function(df = NULL,
                           sizesDf = NULL,
                           cpi = NULL,
                           full = TRUE,
                           ...) {

  #### tests ####
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates
  df$afdm_mg_m2 <- df$n_m2 * df$afdm_mg
  N.ann.list = estimate_ann_stats(df,
                     var = 'n_m2')
  B.ann.list = estimate_ann_stats(df,
                                  var = 'afdm_mg_m2')
 df$afdm_mg_m2 <- NULL

  #### calculate SAMPLE annual production ####
  # Create a matrix with these 4 columns: individual length (mm), mean density for all samples throughout year (number m^-2), individual mass (mg AFDM), and biomass (mg AFDM m^-2) for each size class (rows)
  SF <- matrix(0, length(unique(unlist(df$lengthClass))), 4)

  SF[, c(1, 3)] <- c(
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
  SF.int[, 4] <- SF.int[, 3] * nrow(sizesDf)
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
  # Calculate annual production using the cohort production interval (cpi) given in days for this taxon
  P.ann.samp <- P.uncorr.samp * (365 / cpi)

  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,P.uncorr.samp = P.uncorr.samp,
                B.ann.mean = B.ann.list$afdm_mg_m2_mean,
                B.ann.sd = B.ann.list$afdm_mg_m2_sd,
                N.ann.mean = N.ann.list$n_m2_mean,
                N.ann.sd = N.ann.list$n_m2_sd))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                B.ann.samp = B.ann.list$afdm_mg_m2_mean,
                N.ann.samp = N.ann.list$n_m2_mean))
  }
}
