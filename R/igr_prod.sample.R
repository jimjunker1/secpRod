#' @name igr_prod.sample
#' @title igr_prod.sample
#' @description This function calculates taxa production based on the instantaneous growth method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the
#' @param massValue character string identifying the column name of the mass value
#' @param massLabel character string identifying the column name of the mass value
#' @param full logical. should the full summary be returned with mean and sd.
#' @param ... additional arguments passed to function, including variables to predict growth rate from growth function
#' @return list object with taxa summary of the sampled data
#' @export
igr_prod.sample <- function(df = NULL,
                           sizesDf = NULL,
                           massValue = NULL,
                           massLabel = NULL,
                           full = TRUE,
                           ...) {

  #### tests ####
  #declare globals
  #### GUTS of function ####
  # calculate mean biomass and abundance across all dates
  df[[massLabel]] <- df$n_m2 * df[[massValue]]
  N.ann.list = estimate_ann_stats(df,
                                  var = 'n_m2')
  B.ann.list = estimate_ann_stats(df,
                                  var = massLabel)
  df[[massLabel]] <- NULL

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
