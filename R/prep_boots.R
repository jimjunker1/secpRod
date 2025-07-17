#' @description
#' This function prepares bootstrap samples.
#'
#' @title prep_boots
#' @param df data.frame. A dataframe of the species size, mass, and frequency data.
#' @param bootNum integer. The number of bootstrapped data sets that should be created.
#' @importFrom tidyr pivot_wider
#' @importFrom stats var
#' @export

prep_boots <- function(df = NULL,
                       bootNum = bootNum) {
  ### tests ###
  samp_num <- unlist(aggregate(df$repID, by = list(df$dateID), FUN = function(x) length(unique(x)))$x)
  if(is.na(var(samp_num))){
    warning("Warning: only one sample is present with individuals.")
    } else if (var(samp_num) != 0){
    warning("Warning: sample numbers are not equal among dates")
    }

  # allocate vectors
  bootList <- vector(mode = "list", length = bootNum)

  # create the boot lists
  bootList <- lapply(bootList,
    FUN = function(x) {
      do.call(
        rbind,
        mapply(FUN = function(x,y) {
          x[which(x$repID %in% sample(1:y, y, replace = TRUE)), ]
        }, x = split(df, df$dateID), y = samp_num, SIMPLIFY = FALSE)
      )
    }
  )

  return(bootList)
}
