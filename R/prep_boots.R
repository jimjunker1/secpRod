#' @description
#' This function prepares bootstrap samples.
#'
#' @title prep_boots
#' @param df data.frame. A dataframe of the species size, mass, and frequency data.
#' @param bootNum integer. The number of bootstrapped data sets that should be created.
#' @importFrom stats var
#' @export

prep_boots <- function(df = NULL,
                       bootNum = bootNum,
                       dateCol = 'dateID',
                       repCol = 'repID') {
  ### tests ###
  samp_num <- unlist(aggregate(df[[repCol]], by = list(df[[dateCol]]), FUN = function(x) length(unique(x)))$x)
  if(is.na(var(samp_num))){
    warning("Warning: only one sample is present with individuals.")
    } else if (var(samp_num) != 0){
    warning("Warning: sample numbers are not equal among dates")
    }

  # allocate vectors
  bootList <- vector(mode = "list", length = bootNum)

  # create the boot lists
  ## for all 1:bootNum
  ## sample
  bootList <- lapply(bootList,
    FUN = function(x) {
      do.call(
        # bind all dates together
        rbind,
        # sample replicates within dates
        mapply(FUN = function(x,y,z) {
          ## create within date replicate sample vector of samp_num length
          repSamps = sample(z,y, replace = TRUE)
          ## create a recoded sample identifier to distinguish resampled replicates
          recode = 1:y
          ## sample recode values. keeping original replicate as 1_{old repID}, 2_{old repID}, etc.
          newRepID = as.list(paste(recode, repSamps, sep = "_"))
          ## select the samples into a list of samp_num
          sampList = lapply(repSamps, function(e) x[which(x[[repCol]] %in% e),])
          ## we pass the list of new samples and recoded repIDs
          ## recode the repCol to new repID
          ## bind them together again
          # doing this in a loop because I need it to work and below is running me around
          # for(i in 1:length(sampList)){
          #   sampList[[i]][[repCol]] <- newRepID[[i]]
          # }
          sampList <- Map(function(a,b){
            a[[repCol]] <- b
            a
          }, sampList, newRepID)
          do.call(rbind.data.frame,sampList)
          # mapply(FUN = function(a,b){
          #   a = data.frame(a)
          #   a[[repCol]] <- rep(b, nrow(a))
          #   return(a)
          # },
          # # list of new samples
          # a = sampList,
          # # vector of sample recode values
          # b = newRepID
          # )
        },
        # full taxaSampleMassList split by date
        x = split(df, df[[dateCol]]),
        # # of replicates for each date
        y = as.list(samp_num),
        #list of repIDs for each date
        z = lapply(split(df, df[[dateCol]]), function(x) unique(unlist(x[[repCol]]))),
        SIMPLIFY = FALSE
      ))
    }
  )

  return(bootList)
}
