#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param df
#' @param bootNum
prep_boots <- function(df = taxaSampleListMass,
                       bootNum = bootNum) {

  ### tests ###
  samp_num = unlist(aggregate(df$repID, by = list(df$dateID), FUN = function(x) length(unique(x)))$x)
  if(var(samp_num) !=0 ) warning("Warning: sample numbers are not equal among dates")

  #allocate vectors
  bootList = vector(mode = 'list', length = bootNum)

  # create the boot lists
  bootList = lapply(bootList,
                    FUN = function(x) do.call(rbind,
                                              lapply(split(df, df$dateID),
                                                     function(x) x[sample(nrow(x), unique(samp_num), replace = TRUE), ])))

  return(bootList)

}
