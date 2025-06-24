#' @description
#' This function calculates production using the instantaneous growth rate method.
#' @details
#' Additional details...
#' @title calc_prod_igr
#' @param taxaSampleListMass description
#' @param taxaInfo The taxa info data.frame
#' @param bootNum integer. The number of bootstrap samples to produce.
#' @param dateDf data frame of date information with external predictors for each month. There should be a column name identical to all variables in the growth equation found in taxaInfo data.frame.
#' @param taxaSummary string of \code{'short'}, \code{'full'}, or \code{'none'}. What type of summary information should be returned.
#' @param wrap logical. Should the dates wrap to create a full year?
#' @param massValue string. What is the mass value and units of the production
#' @param massLabel string. What label should the output units be. It is possible this will default to 'massValue' in the future.
#' @param bootList list. This is the bootstrapped samples passed from `calc_production()`
#' @param envData an optional data frame (or coercible) of additional variables that are called to feed into growth equations to estimate growth rates. The variable column names must match exactly the terms in the growth equation from taxaInfo.
#' @param ... additional arguments to pass to the function
#' @returns returns a list of 2 objects
#' @importFrom dplyr count
#' @import stats
#' @export

calc_prod_igr <- function(taxaSampleListMass= NULL,
                          taxaInfo = NULL,
                          bootNum = NULL,
                          dateDf = NULL,
                          taxaSummary = 'full',
                          wrap = TRUE,
                          massValue = NULL,
                          massLabel = NULL,
                          bootList = NULL,
                          envData = NULL,
                          ...) {

  ## tests ##
  # check for environmental data

  # check that the dates match

  # check for parsed variables and their presence in envData
  # - parse all the unique variables
  ## input the growth formula
  growthFormula <- as.formula(paste0(gsub(" ","",taxaSubInfo$growthForm)))
  # - confirm if g_d is present, check for other contributed variables e.g., massValue, lengthClass,etc.
  ##
  growthUnits <- formula.tools::lhs(growthFormula)
  if(!grepl("g_d", growthUnits)){stop("Error: Currently, growth formula must contain `g_d` as growth rate variable name.")}

  # - check for novel variables, e.g., temperature, OM, etc. and their presence in envData
  # convert the length-to-mass formula to a formula
  massFormula <- as.formula(paste0(gsub(" ", "", taxaSubInfo$massForm)))
  # if(!as.formula(massFormula)) stop("The length-to-mass formula, 'massForm', cannot be coercied into a formula class.")
  # get the mass units from the LHS of massFormula
  massUnits <- formula.tools::lhs(massFormula)
  # get the variables from the RHS of massFormula
  massRHS <- formula.tools::rhs(massFormula)
  # detect if polynomial expression exists
  caret_present <- any(grepl("\\^", massRHS, ignore.case = TRUE))
  # if (caret_present) warning("R cannot parse functions with '^'.)
  if (!any(grepl("lengthClass", massRHS, ignore.case = FALSE))) stop("`lengthClass` is not present in the length-mass (L-M) formula. This variable must be present to convert from length to mass with user-defined L-M function.")
  ## The code below is a work in progress to parse the formula when `^` is present. This is a feature for later. I will likely remove this note and code and add it to a development branch in the near future.
  # if(caret_present){
  #   # deparse to convert 'call' class to character
  #   charRHS = deparse(massRHS)
  #   # split on punctuation to extract just variable names
  #   allVars = sapply(unlist(strsplit(charRHS, "[[:punct:]]")), trimws)
  #   # remove any blank spaces from the vector
  #   allVars = otherVars[grepl( "\\w", otherVars)]
  #   # detect how many variables do not have columns
  #   if(sum(massRHS %ni% colnames(taxaSubInfo)) > 1) stop("There are multiple undefined terms in the length-to-mass equation. There should only be a one (1) corresponding to lengthClass.")
  #   lengthVar = allVars[allVars %ni% colnames(taxaSubInfo)]
  #   names(lengthVar) <- 'lengthClass'
  #
  #   powerVars = unlist(strsplit(gsub("^.+(\\w\\^\\w).+$","\\1", charRHS),"\\^"))
  #   collapsePowerVars = paste(powerVars, collapse = "")
  #   removeCaret = gsub("\\^"," ", charRHS)
  #   replaceCaret = gsub(colPowerVars, powerVars, removeCaret)
  #
  #
  #   varRHS = rhs.vars(newRHS)
  # }
  # confirm all variables are present in the taxaInfo sheet
  # deparse to convert 'call' class to character
  charRHS <- deparse(massRHS)
  # split on punctuation to extract just variable names
  allVars <- sapply(unlist(strsplit(charRHS, "[[:punct:]]")), trimws)
  # remove any blank spaces from the vector
  allVars <- allVars[grepl("\\w", allVars)]
  if (!any(grepl("lengthClass", allVars, ignore.case = FALSE))) stop("lengthClass must be a term in the length-to-mass equation.")
  # remove "poly" for polynomials
  allVars <- allVars[allVars %ni% "lengthClass"]
  # detect how many variables do not have columns
  if (sum(allVars %ni% names(taxaSubInfo)) > 1) stop("There are multiple undefined terms in the length-to-mass equation. There should only be a one (1) corresponding to lengthClass.")

  # fill in the named variables L-M equation with their taxon-specific values
  ## create a named list of all variables
  allVars <- setNames(allVars, nm = taxaSubInfo[, names(taxaSubInfo) %in% allVars])
  ## parse and replace the RHS of the L-M function call
  newRHS = sapply(charRHS, function(x) replace_formula_terms(x, termList = allVars))
  ## Recombine the LHS & RHS to a new function
  newMassFormula = as.formula(paste0(massUnits,"~",newRHS))
  massUnits = lhs(newMassFormula)
  massRHS = rhs(newMassFormula)

  ## end tests ##
  speciesName = unique(taxaSampleListMass$taxonID)
  # ## function prep ##
  # ### make a list of key variables to pass to sample function
  funcList = list(
    df = taxaSampleListMass,
    # sizesDf = unique(taxaSampleListMass[, c("lengthClass", rev(names(taxaSampleListMass))[1])])
    sizesDf = unique(taxaSampleListMass[, c("lengthClass", eval(massValue))]),
    massValue = massValue,
    massLabel = massLabel,
    envData = envData
  )

  # calculate the production from the full sample
  if(is.numeric(taxaInfo$growthForm)){
    if(length(taxaInfo$growthForm) == 1){
      growthForm = unlist(taxaInfo$growthForm)
    }
  }

  taxaCPI <- mean(c(taxaInfo$min.cpi, taxaInfo$max.cpi))
  funcList = c(funcList, list(cpi = taxaCPI))
  P.samp = do.call(igr_prod.sample, args = funcList)

  # prep boots
  # bootList = prep_boots(df = taxaSampleListMass,
  #                        bootNum = bootNum)

  P.boots = mapply(FUN = igr_prod.sample,
                   df = bootList,
                   sizesDf = lapply(1:bootNum, function(x) funcList$sizesDf),
                   massValue = massValue,
                   massLabel = massLabel,
                   full = FALSE)

  #### create SAMPLE information to export as summary ####
  # summarise sample sizes across dates
  sampDatesInfo <- setNames(unique(aggregate(taxaSampleListMass[c("repID")], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), count)[c(1, 3)]), c("dateID", "N"))
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[["N"]] <- NA
    sampDatesInfo <- rbind(sampDatesInfo, temp)
  }
  # summarise the sample abundance N across all dates and size classes
  Nmean <- setNames(cleanAggDf(aggregate(taxaSampleListMass, by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", "n_m2_mean"))
  Nsd <- setNames(cleanAggDf(aggregate(taxaSampleListMass["n_m2"], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), sd, na.rm = TRUE)), nm = "n_m2_sd")
  Nbind <- cbind(Nmean, Nsd)
  Nbind$lengthClass <- factor(Nbind$lengthClass, levels = unique(Nbind$lengthClass))
  NmeanTab <- as.data.frame.matrix(xtabs(n_m2_mean ~ dateID + lengthClass, Nbind))
  NsdTab <- as.data.frame.matrix(xtabs(n_m2_sd ~ dateID + lengthClass, Nbind))
  NdatesInfo <- setNames(aggregate(Nmean["n_m2_mean"], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", "n_m2_mean"))
  # if wrap equals true create another
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[["n_m2_mean"]] <- mean(c(NdatesInfo[1, "n_m2_mean"], NdatesInfo[nrow(NdatesInfo), "n_m2_mean"]))
    NdatesInfo <- rbind(NdatesInfo, temp)
  }
  # summarise the sample biomasses across all dates and size classes
  # create the sizeclass biomass for all rows
  taxaSampleListMass[[massLabel]] <- unlist(taxaSampleListMass[, "n_m2"]) * unlist(taxaSampleListMass[, massValue])
  # do the aggregating
  Bmean <- setNames(cleanAggDf(aggregate(taxaSampleListMass[c("dateID", "lengthClass", massLabel)], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), mean, na.rm = TRUE)), nm = c("dateID", "lengthClass", paste0(massLabel, "_mean")))
  Bsd <- setNames(cleanAggDf(aggregate(taxaSampleListMass[massLabel], by = list(taxaSampleListMass$dateID, taxaSampleListMass$lengthClass), sd, na.rm = TRUE)), nm = paste0(massLabel, "_sd"))
  Bbind <- cbind(Bmean, Bsd)
  Bbind$lengthClass <- factor(Bbind$lengthClass, levels = unique(Bbind$lengthClass))
  meanBform <- as.formula(paste0(massLabel, "_mean ~ dateID + lengthClass"))
  sdBform <- as.formula(paste0(massLabel, "_sd ~ dateID + lengthClass"))
  BmeanTab <- as.data.frame.matrix(xtabs(meanBform, Bbind))
  BsdTab <- as.data.frame.matrix(xtabs(sdBform, Bbind))
  BdatesInfo <- setNames(aggregate(Bmean[paste0(massLabel, "_mean")], by = list(Nmean$dateID), sum, na.rm = TRUE), nm = c("dateID", paste0(massLabel, "_mean")))
  # if wrap equals true
  if (wrap) {
    temp <- data.frame(dateID = dateDf[nrow(dateDf), "dateID"])
    temp[[eval(paste0(massLabel, "_mean"))]] <- mean(c(BdatesInfo[1, eval(paste0(massLabel, "_mean"))], BdatesInfo[nrow(BdatesInfo), eval(paste0(massLabel, "_mean"))]))

    BdatesInfo <- rbind(BdatesInfo, temp)
  }
  # create the full summary
  datesInfo <- Reduce(function(x, y) merge(x, y, all = TRUE), list(sampDatesInfo, NdatesInfo, BdatesInfo))

  if(taxaSummary == "none"){

  } else if (taxaSummary == "full") {
    # # create a list for output
    taxaSummary <- list(
      summaryType = "full",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.samp$P.ann.samp,
      P.uncorr.samp = P.samp$P.uncorr.samp,
      cpi = taxaCPI,
      meanN = mean(unlist(datesInfo$n_m2_mean)),
      meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
      meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
      Nmean = NmeanTab,
      Nsd = NsdTab,
      Bmean = BmeanTab,
      Bsd = BsdTab,
      datesInfo = datesInfo
    )
  } else if(taxaSummary == "short"){
    taxaSummary <- list(
      summaryType = "short",
      taxonID = taxaInfo$taxonID,
      method = "sf",
      P.ann.samp = P.samp$P.ann.samp,
      cpi = taxaCPI,
      meanN = mean(unlist(datesInfo$n_m2_mean)),
      meanB = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])),
      meanIndMass = mean(unlist(datesInfo[[eval(paste0(massLabel, "_mean"))]])) / mean(unlist(datesInfo$n_m2_mean)),
      datesInfo = datesInfo
    )
  }
  #   assign(taxaInfo$taxonID, list())
  # # add taxainformation
  #   assign(taxaInfo$taxonID,
  #          within(eval(as.symbol(taxaInfo$taxonID)),{
  #            taxaInfo <- taxaInfo
  #            }
  #            )
  #          )

  return(assign(speciesName, list(P.boots = P.boots,
                                  taxaSummary = taxaSummary)))

}
