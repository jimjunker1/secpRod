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
#' @param minGrowth numeric. a minimum growth rate value for when growth is negative.
#' @param wrap logical. should the calculations be wrapped by adding an additional date to make a full year?
#' @param full logical. should the full summary be returned with mean and sd.
#' @param ... additional arguments passed to function, including variables to predict growth rate from growth function
#' @return list object with annual production, mean biomass, and mean abundance of the sample
#' @export
igr_prod.sample <- function(df = NULL,
                           dateDf = NULL,
                           growthDf = NULL,
                           lengthValue = NULL,
                           massValue = "mass",
                           abunValue = "density",
                           dateCol = "dateID",
                           repCol = "repID",
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


  #### calculate SAMPLE annual production ####
  # merge the df with the growthDf and calculate size-specific production
  df <- merge(df, growthDf)#, by = c(dateCol, massValue,lengthValue))
  if(wrap == TRUE){
    # get the detilas of first and final dates with samples
    firstLast <- dateDf[c(1,(nrow(dateDf)-1)),]
    wrapSamples <- dplyr::filter(df, df[[dateCol]] %in% unique(unlist(firstLast[[dateCol]])))
    wrapSamples$taxonID <- NULL;wrapSamples[[repCol]] <- NULL
    wrapForm_1 <- paste0(".~", paste0(c(dateCol,lengthValue,massValue), collapse = "+"))
    wrapAgg_first <- stats::aggregate(formula(wrapForm_1), data = wrapSamples, na.rm = TRUE,  FUN = mean)
    wrapAgg_first[[dateCol]] <- NULL
    wrapAgg_first[['int_days']] <- 1
    wrapForm_2 <- gsub(paste0(dateCol,"+."), "",wrapForm_1)
    wrapProduction <- stats::aggregate(formula(wrapForm_2), data = wrapAgg_first, na.rm = TRUE, FUN = mean)
    wrapProduction$production <- wrapProduction[["biomass"]] * wrapProduction[['g_d']]
    wrapProduction[[dateCol]] <- dateDf[nrow(dateDf), dateCol]
    wrapProduction <- stats::aggregate(formula(paste0("production~",dateCol)), data = wrapProduction, FUN = sum, na.rm = TRUE)
  }
  df$production <- df[["biomass"]] * df[["g_d"]] * df[["int_days"]]
  dfRepAgg <- stats::aggregate(formula(paste0("production~",dateCol,"+",repCol)), data = df, sum, na.rm = TRUE)
  dfDateAgg <- stats::aggregate(formula(paste0("production~",dateCol)), data = dfRepAgg, mean)
  if(wrap == TRUE){
    dfDateAgg <- rbind(dfDateAgg, wrapProduction)
  }
  P.ann.samp <- sum(dfDateAgg$production)


  if(full == TRUE){
    return(list(P.ann.samp = P.ann.samp,
                B.ann.mean = B.ann.list[["biomass_mean"]],
                B.ann.sd = B.ann.list[["biomass_sd"]],
                N.ann.mean = N.ann.list[[paste0(abunValue,"_mean")]],
                N.ann.sd = N.ann.list[[paste0(abunValue,"_sd")]],
                dfDateAgg = dfDateAgg))
  } else{
    return(list(P.ann.samp = P.ann.samp,
                B.ann.samp = B.ann.list[["biomass_mean"]],
                N.ann.samp = N.ann.list[[paste0(abunValue,"_mean")]]))
  }
}
