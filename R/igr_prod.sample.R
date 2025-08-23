#' @name igr_prod.sample
#' @title igr_prod.sample
#' @description This function calculates taxa production based on the instantaneous growth method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the sizes
#' @param lengthValue string of the column name containing the length class measurements
#' @param massValue character string identifying the column name of the mass value
#' @param abunValue character string identifying the column name of the density value
#' @param dateCol character string identifying the column name of the sample date information. This is distinguished from *Value parameters in that Values may be used to maintain units provenance in the future. This may also change.
#' @param repCol character string identifying the column name of the replicate information.
#' @param wrap logical. should the calculations be wrapped by adding an additional date to make a full year?
#' @param full logical. should the full summary be returned with mean and sd.
#' @param ... additional arguments passed to function, including variables to predict growth rate from growth function
#' @return list object with annual production, mean biomass, and mean abundance of the sample
#' @export
igr_prod.sample <- function(df = NULL,
                           dateDf = NULL,
                           lengthValue = NULL,
                           massValue = massValue,
                           abunValue = abunValue,
                           dateCol = dateCol,
                           repCol = repCol,
                           wrap = FALSE,
                           full = TRUE,
                           ...) {

  #### tests ####
  #declare globals
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates
  df[["biomass"]] <- df[[abunValue]] * df[[massValue]]
  N.ann.list = estimate_ann_stats(df, var = abunValue,
                                  massValue = massValue,
                                  abunValue = abunValue,
                                  dateCol = dateCol,
                                  repCol = repCol,
                                  wrap = wrap)

  B.ann.list = estimate_ann_stats(df, var = "biomass",
                                  massValue = massValue,
                                  abunValue = abunValue,
                                  dateCol = dateCol,
                                  repCol = repCol,
                                  wrap = wrap)

  # parse the growth equation to predict growth rates

  #### calculate SAMPLE annual production ####
  # Create a matrix with these 4 columns: individual length (mm), mean density for all samples throughout year (number m^-2), individual mass (mg AFDM), and biomass (mg AFDM m^-2) for each size class (rows)

  return(NULL)
  # if(full == TRUE){
  #   return(list(P.ann.samp = P.ann.samp,P.uncorr.samp = P.uncorr.samp,
  #               B.ann.mean = B.ann.list$afdm_mg_m2_mean,
  #               B.ann.sd = B.ann.list$afdm_mg_m2_sd,
  #               N.ann.mean = N.ann.list$n_m2_mean,
  #               N.ann.sd = N.ann.list$n_m2_sd))
  # } else{
  #   return(list(P.ann.samp = P.ann.samp,
  #               B.ann.samp = B.ann.list$afdm_mg_m2_mean,
  #               N.ann.samp = N.ann.list$n_m2_mean))
  # }
}
