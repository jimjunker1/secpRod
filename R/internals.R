#' @description
#' This function is the opposite of \code{\%in\%} in that it finds items not in a vector
#' @title \%ni\%
#' @param x vector or NULL: the values to exclude. Long vectors are supported.
#' @param table vector or NULL: the values to be excluded against. Long vectors are not supported.
#' #' @returns A logical vector, indicating if a match was not located for each element of x: thus the values are TRUE or FALSE and never NA.
`%ni%` <- function(x, table) {!(x %in% table)}
# '%ni%' <- Negate('%in%')


#' @title inv.logit
#' @param eta from mixtools
inv.logit <- binomial()$linkinv

#'
#'
factorise <- function(x) {
  x <- length(x)
  if(x == 1){return(1)}
  todivide <- seq(from = 2, to = x)
  out <- todivide[x %% todivide == 0L]
  return(out)
}

#' @title findreps
#' @description This function finds and sets breaks on repeated runs of similar character types. It is used to parse the formula structure and variables in `massForm` for converting lengths to mass in \code{convert_length_to_mass()}.
#' @param x character string of the formula structure
#' @param counter stuff
#' @source https://stackoverflow.com/questions/33155662/find-and-break-on-repeated-runs
findreps <- function(x, counter = NULL){
  if(is.null(counter)){
    counter <- c()
    maxcounter <- 0
  } else {
    maxcounter <- max(counter)
  }
  holding <- lapply(1:length(x), function(y){x[1:y]})
  factors <- lapply(holding, factorise)
  repeats <- sapply(1:length(factors), function(index) {any(sapply(1:length(factors[[index]]), function(zz) {all((rep(holding[[index]][1:(length(holding[[index]])/factors[[index]][zz])], factors[[index]][zz]))==holding[[index]])}))})
  holding <- holding[max(which(repeats))][[1]]
  if(length(holding) == length(x)){
    return(c(counter, rep(maxcounter + 1, length(x))))
  } else {
    counter <- c(counter, rep(maxcounter + 1, length(holding)))
    return(findreps(x[(length(holding) + 1):length(x)], counter))
  }
}


#'
#'
#'
replace_formula_terms <- function(x, termList = NULL, ...) {
  "%ni%" <- Negate("%in%")

  # termList = setNames(names(df), nm = unlist(df[1,]))

  otherVars = unlist(gregexpr("[[:punct:]]", x))
  names(otherVars) = sapply(otherVars, function(a) substring(x,a,a))
  lengthVar = unlist(gregexpr("lengthClass", x))
  names(lengthVar) = "lengthClass"
  termVars = unlist(gregexpr("\\w", x))
  termVars = termVars[termVars %ni% lengthVar:(lengthVar+11)]
  termReps = findreps(diff(termVars))
  termRepsTable = table(termReps)
  ones = as.numeric(names(termRepsTable)[termRepsTable == 1])
  onesPositions = termVars[sapply(ones, function(x) grep(x,termReps))]
  names(onesPositions) = sapply(substring(x, onesPositions, onesPositions), function(a) names(termList)[termList %in% a])
  runs = as.numeric(names(termRepsTable)[termRepsTable > 1])
  if(length(sapply(runs, function(a) grep(a,
                                          termReps))) == 0){
    runsPositions = NA
    onesPositions = termVars
    names(onesPositions) = sapply(substring(x, onesPositions, onesPositions), function(a) names(termList)[termList %in% a])
  } else{
  runsPositions = termVars[sapply(runs, function(a) grep(a, termReps))]}
  runsStarts = min(runsPositions)
  names(runsStarts) = sapply(substring(x,min(runsPositions),(max(runsPositions)+1)), function(a) names(termList)[termList %in% a])

  newRHSForm = paste(names(sort(c(otherVars, lengthVar, onesPositions, runsStarts))), collapse = "")

}

#'
#'
#'
MissLt <- function(x, ratio = 0.5){
  sum(is.na(x))/length(x) <= ratio
}

#'
#'
#'
dateCoercible  <- function(x, addformat = NULL, exactformat = NULL){
  if (is.null(exactformat)){
    format = c("%m/%d/%Y", "%m-%d-%Y","%Y/%m/%d" ,"%Y-%m-%d", addformat)
    y <- as.Date(as.character(x),format= format)
    MissLt(y,ratio = 1-(1/length(y)))}
  else{
    y <- as.Date(as.character(x),format= exactformat)
    MissLt(y,ratio = 1-(1/length(y)))}
}


#'
#'
#'
cleanAggDf = function(df,...){
  good_cols = !grepl("Group\\.\\d{1}|repID", names(df), ignore.case = TRUE)
  filter_cols = Filter(function(x)!all(is.na(x)), df[good_cols])
  if(dim(filter_cols)[2] < 1){
    data.frame(x = NA_real_)
  } else{
    filter_cols
  }
}

#'
#'
#'
estimate_ann_stats = function(df = NULL,
                              var = NULL,
                              wrap = TRUE,
                              massValue = 'afdm_mg',
                              abunValue = 'density',
                              dateCol = 'dateID',
                              repCol = 'repID',
                              ...){
  if(nrow(df) == 0){
    varMeanName = paste0(var,"_mean")
    varSDName = paste0(var,"_sd")
    x = list()
    x[[varMeanName]] <- 0
    x[[varSDName]] <- NA_real_
    return(x)
  }
  # sum across lengthClass
  varDateSum <- stats::setNames(stats::aggregate(df[var], by = list(df[[dateCol]], df[[repCol]]), sum, na.rm = TRUE), nm = c("dateID", "repID", paste0(var,"_sum")))
  # take means of all dates across reps
  varDateSumMean <- stats::setNames(stats::aggregate(varDateSum[grep(var, names(varDateSum))], by = list(varDateSum$dateID), mean, na.rm = TRUE), nm = c("dateID", paste0(var,"_mean")))
  # take sd of all dates across reps
  varDateSumSD <- stats::setNames(stats::aggregate(varDateSum[grep(var, names(varDateSum))], by = list(varDateSum$dateID), sd, na.rm = TRUE), nm = c("dateID", paste0(var,"_sd")))
  # wrap the date around
  if(wrap){
    varSumMeanWrap = (varDateSumMean[1,grep(var, names(varDateSumMean))]+varDateSumMean[nrow(varDateSumMean),grep(var, names(varDateSumMean))])/2
    varSumMeanWrapSD = sqrt(varDateSumSD[1,grep(var, names(varDateSumSD))]^2+varDateSumSD[nrow(varDateSumSD),grep(var, names(varDateSumSD))]^2)
    varDateWrap = as.Date(varDateSumMean[1,dateCol]+364)
    wrapDf = data.frame(varDateWrap)
    names(wrapDf)[1] <- dateCol
    wrapDf[paste0(var,"_mean")] <- varSumMeanWrap
    wrapDf[paste0(var,"_sd")] <- varSumMeanWrapSD
    varDateSumMean <- rbind(varDateSumMean, wrapDf[,c(1,2)])
    varDateSumSD <- rbind(varDateSumSD, wrapDf[,c(1,3)])
  }

  # estimate the annual mean
  varMean = mean(unlist(varDateSumMean[grep(var, names(varDateSumMean))]), na.rm = TRUE)
  # estimate the SD of the annual mean
  varSD = sqrt(sum(unlist(varDateSumSD[grep(var, names(varDateSumSD))])^2, na.rm = TRUE))

  varMeanName = paste0(var,"_mean")
  varSDName = paste0(var,"_sd")
  x = list()
  x[[varMeanName]] <- varMean
  x[[varSDName]] <- varSD
  return(x)
}

#'
#'
vers_count = function (..., condition = (function(x) TRUE))
{
  data <- c(...)
  result <- sum(sapply(data, function(x) if (condition(x)) 1 else 0))
  return(result)
}

#'
#'
string_agg = function(df, colString = NULL){
  if(length(colString) > 1){
    firstString = paste0("cbind(",paste0(colString, collapse = ","),")~")
    secondString = paste0('taxonID+dateID+repID+',eval(massValue))
  } else{
    firstString = paste0(colString,"~")
    secondString = paste0('taxonID+dateID+repID+',eval(massValue))
  }
  aggregate(formula(paste0(firstString, secondString)), data = df, FUN = vers_count)
}

#'
#'
#'
get_bin_widths = function(x){
  boundsString = sapply(unlist(x), FUN = function(x) gsub("\\[|\\]|\\(","",x))
  boundsSplit = stats::setNames(data.frame(t(sapply(sapply(boundsString, function(x) strsplit(x,",")), as.numeric))), nm = c('bin_min','bin_max'))
  boundsSplit$midpoint = as.character(round(apply(boundsSplit,1, FUN = function(x) mean(c(x[1],x[2]))),3))
  rownames(boundsSplit) <- NULL
  return(boundsSplit)
}

#'
# recode_boot_dates = function(a = NULL, b = NULL){
#
# }
#'
#'
# calc_Linf = function(df,...){
#   dataList = make_standata(log(count)~lengthClass, data = df)
#   modelCode = make_stancode(log(count)~lengthClass,
#                             family = gaussian(),
#                             prior = set_prior('normal(0,1)', class = 'b', ub = 0),
#                             data = df)
#
#   l_nLM = stan(model_code = modelCode, data = dataList,
#                chains = 4, iter = 2000, warmup = 1000,  thin = 1)
#
#   int = unlist(extract(l_nLM)['Intercept'])
#   b = unlist(extract(l_nLM)['b'])
#   Linf = mapply(FUN = function(a,b) a/(1-b), int, b)
#   ZK = sapply(b, function(b) b/(1-b))
#
#   return(list(model = l_nLM,
#               Linf = Linf,
#               ZK = ZK))
# }
