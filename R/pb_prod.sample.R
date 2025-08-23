#' @name pb_prod.sample
#' @title pb_prod.sample
#' @description This function calculates sample production based on the production:biomass method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param massValue character string identifying the column name of the mass value
#' @param abunValue character string identifying the column name of the abundance value
#' @param dateCol character string identifying the column name of the date column
#' @param repCol character string identifying the column name of the replicate column
#' @param pb numeric. The production:biomass ratio.
#' @param wrap logical. should the calculations be wrapped by adding an additional date to make a full year?
#' @param full logical. should the full summary be returned with mean and sd
#' @param ... additional arguments passed to function
#' @return list object with annual production, mean biomass, and mean abundance of the sample
#' @export
pb_prod.sample <- function(df = NULL,
                           massValue = 'mass',
                           abunValue = 'density',
                           dateCol = 'dateID',
                           repCol = 'repID',
                           pb = NULL,
                           wrap = FALSE,
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
 # df[['biomass']] <- NULL
 if(B.ann.list[["biomass_mean"]] == 0){
   if(full == TRUE){
     return(list(P.ann.samp = 0,
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
  # Calculate annual production using the production to biomass ratio given for this taxon
  P.ann.samp <- B.ann.list[["biomass_mean"]]*pb
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
