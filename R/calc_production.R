#' @name calc_production
#' @description This is the main function of the secpRod package. It will calculate all secondary production for all groups based on the methods described in the taxa information object. Depending on input values varying summaries
#' @param sizeInfo dataframe (or coercible); The size abundance data for each species and sampling date
#' @param infoCols integer vector; Any columns in the sizeInfo object that are not the taxonomic ID, sampling data, or size class columns
#' @param taxInfo dataframe (or coercible); The taxonomic information for calculating secondary production. This must include a taxonomic ID column with the same name as that of `sizeInfo`, other information
#'
