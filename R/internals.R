#'
#'
#'
factorise <- function(x) {
  x <- length(x)
  if(x == 1){return(1)}
  todivide <- seq(from = 2, to = x)
  out <- todivide[x %% todivide == 0L]
  return(out)
}

#'
#'
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
  runsPositions = termVars[sapply(runs, function(a) grep(a, termReps))]
  runsStarts = min(runsPositions)
  names(runsStarts) = sapply(substring(x,min(runsPositions),(max(runsPositions)+1)), function(a) names(termList)[termList %in% a])

  newRHSForm = paste(names(sort(c(otherVars, lengthVar, onesPositions, runsStarts))), collapse = "")

}
