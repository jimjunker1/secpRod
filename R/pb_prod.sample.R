#' @name pb_prod.sample
#' @title pb_prod.sample
#' @description This function calculates taxa production based on the production:biomass method
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param sizesDf a data.frame of the
#' @param massValue character string identifying the column name of the mass value
#' @param massLabel character string identifying the column name of the mass value
#' @param full logical. should the full summary be returned with mean and sd
#' @param pb numeric. The production:biomass ratio.
#' @param ... additional arguments passed to function
#' @return list object with taxa summary of the sampled data
#' @export
pb_prod.sample <- function(df = NULL,
                           sizesDf = NULL,
                           massValue = NULL,
                           massLabel = NULL,
                           pb = NULL,
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
  P.ann.samp <- B.ann.list[[paste0(massLabel,"_mean")]]*pb
}
  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,
                B.ann.mean = B.ann.list[[paste0(massLabel,"_mean")]],
                B.ann.sd = B.ann.list[[paste0(massLabel,"_sd")]],
                N.ann.mean = N.ann.list$n_m2_mean,
                N.ann.sd = N.ann.list$n_m2_sd))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                B.ann.samp = B.ann.list[[paste0(massLabel,"_mean")]],
                N.ann.samp = N.ann.list$n_m2_mean))
  }

}
