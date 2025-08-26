#' @title convert_length_to_mass
#' @description A function to convert species-specific lengths to mass based on a user-provided length-mass equation form and variable values
#' @param taxaSampleList data.frame in long format for a single taxa. The data.frame should contain a species identifier column `taxonID` and a column of length bin categories lengthValue.
#' @param taxaInfo a data.frame of the information to convert length to mass for all taxa. The taxa specified in `taxaSampleList` will be subset from here. This data.frame must contain a `taxonID` column, the length-to-mass equation formula, `massForm`, which must contain `lengthClass` as a variable (e.g., `afdm_mg~a*lengthClass^b`). Additional columns are necessary based on the length-mass formula. All other non-`lengthClass` variables on the right hand side (RHS) must have unique columns named the variable name. For example, the above formula structure, `afdm_mg~a*lengthClass^b`, the RHS is `a*lengthClass^b`. `lengthClass` is a required column, but optionally necessary columns are `a` and `b` for the other variables. The formula will be parsed and species-specific `a` and `b` coefficients will be inserted for conversion.
#' @param lengthValue a character string identifying the column name for the length class. The data should be numeric or coercible
#' @param massValue a character string to name the mass variable. This should not contain any punctuation (e.g., .,_,-, etc.).
#' @param reduce logical. If TRUE (default) the mass column will be added to `taxaSampleList`. The name of the mass column will be parsed from the left hand side (LHS) of the `massForm` provided in `taxaInfo`. For example, in `mass~a*lengthValue^b` the mass column will be named based on the massValue column. This data.frame is returned if reduce == TRUE
#' @param ... additional arguments passed to function
#' @returns taxaSampleList with the mass column added.
#' @importFrom formula.tools rhs lhs
#' @importFrom rlang :=
#' @importFrom stats as.formula
#' @importFrom stats setNames
#' @export
convert_length_to_mass <- function(taxaSampleList = NULL, taxaInfo = NULL, lengthValue = NULL, massValue = NULL, reduce = TRUE, ...) {
  if (is.null(taxaSampleList)) stop("Error: No sample information provided.")
  if (is.null(taxaInfo)) stop("Error: No taxonomic information provided.")
  if (any(is.na(suppressWarnings(as.numeric(unique(taxaSampleList$lengthClass)))))) stop("Error: Non-numeric values in `lengthClass`. Must be numeric or coercible.")
  if (is.na(as.character(massValue))) stop("Error: the `massValue` variable must be a character or coercible")
  if (grepl("[[:punct:]]", massValue)) stop("Error: `massValue` must not contain punctuation. This will break parsing of growth formula in the IGR method.")

  taxonID <- unique(taxaSampleList$taxonID)
  if (length(taxonID) > 1) warning("length(taxonID) > 1")
  # subset the taxaInfo object to a single taxa
  taxaSubInfo <- taxaInfo[which(taxaInfo$taxonID == taxonID), ]
  if (any(is.null(taxaSubInfo))) stop(paste("Error: No taxonomic information available for", taxonID, ". Check for correct spelling in sampleInfo and taxaInfo."))

  # convert the length-to-mass formula to a formula
  massFormula <- stats::as.formula(paste0(gsub(" ", "", taxaSubInfo$massForm)))
  # if(!stats::as.formula(massFormula)) stop("The length-to-mass formula, 'massForm', cannot be coercied into a formula class.")
  # get the mass units from the LHS of massFormula
  massUnits <- formula.tools::lhs(massFormula)
  # get the variables from the RHS of massFormula
  massRHS <- formula.tools::rhs(massFormula)
  # detect if polynomial expression exists
  caret_present <- any(grepl("\\^", massRHS, ignore.case = TRUE))
  # if (caret_present) warning("R cannot parse functions with '^'.)
  if (!any(grepl(lengthValue, massRHS, ignore.case = FALSE))) stop(paste0(lengthValue," is not present in the length-mass (L-M) formula. This variable must be present to convert from length to mass with user-defined L-M function."))
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
  if (!any(grepl(lengthValue, allVars, ignore.case = FALSE))) stop("lengthClass must be a term in the length-to-mass equation.")
  # remove "poly" for polynomials
  allVars <- allVars[allVars %ni% lengthValue]
  # detect how many variables do not have columns
  if (sum(allVars %ni% names(taxaSubInfo)) > 1) stop("There are multiple undefined terms in the length-to-mass equation. There should only be a one (1) corresponding to lengthValue.")

  # fill in the named variables L-M equation with their taxon-specific values
  ## create a named list of all variables
  allVars <- stats::setNames(allVars, nm = taxaSubInfo[, names(taxaSubInfo) %in% allVars])
  ## parse and replace the RHS of the L-M function call
  newRHS = sapply(charRHS, function(x) replace_formula_terms(x, termList = allVars))
  ## Recombine the LHS & RHS to a new function
  newMassFormula = stats::as.formula(paste0(massUnits,"~",newRHS))
  massUnits = lhs(newMassFormula)
  massRHS = rhs(newMassFormula)

  # Mutate length to mass from parsed L-M equation
  if (reduce) {
    return(taxaSampleList %>%
             dplyr::mutate(!!massUnits := !!massRHS))
  }
}
