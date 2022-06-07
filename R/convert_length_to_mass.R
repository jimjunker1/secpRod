#'
#'
#'
convert_length_to_mass <- function(taxaSampleList = NULL, infoCols = NULL, taxaInfo = NULL, reduce = TRUE, ...) {
  library(formula.tools)
  if (is.null(taxaSampleList)) stop("No sample information provided.")
  if (is.null(taxaInfo)) stop("No taxonomic information provided.")

  taxonID <- unique(taxaSampleList$taxonID)
  if (length(taxonID) > 1) warning("length(taxonID) > 1")
  # subset the taxaInfo object to a single taxa
  taxaSubInfo <- taxaInfo[which(taxaInfo$taxonID == taxonID), ]
  if (any(is.null(taxaSubInfo))) stop(paste("No taxonomic information available for", taxonID, ". Check for correct spelling in sampleInfo and taxaInfo."))

  # convert the length-to-mass formula to a formula
  massFormula <- as.formula(gsub(" ", "", taxaSubInfo$massForm))
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

  # Mutate length to mass from parsed L-M equation
  if (reduce) {
    return(taxaSampleList %>%
             dplyr::mutate(!!massUnits := !!massRHS))
  }
}
