#' @description
#' This function creates an additional date in the data set to make a full annual cycle
#' @title wrap_dates
#' @param df a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param envData (optional) a data.frame of the environmental data. The number of rows should be the same as the number of sampling dates.
#' @param dateCol character string identifying the column name of the sampling dates
#' @param wrapDate logical. should we wrap the dates to make full year?
#' @return data.frame with sampling dates and intervals, or number of days between
#' @export
wrap_dates <- function(df = NULL, envData = NULL, dateCol = NULL, wrapDate = TRUE) {
  if(is.null(dateCol)){
    # check if any columns are date
    if(any(unlist(lapply(df, function(x) inherits(x, c("Date","POSIXt")))))){
      # extract the date column, unlist to vector and maintain class
      dateCol = names(df[unlist(lapply(df, function(x) inherits(x, c('Date', 'POSIXt'))))])
      date1 = unname(do.call("c",unique(df[unlist(lapply(df, function(x) inherits(x, c('Date', 'POSIXt'))))])))
      # convert to julian dates. Origin is default "1970-01-01"
      julian1 = julian(date1, origin = as.Date("1970-01-01"))

      t = sort(union(julian1, NULL))

      #Create a vector of durations (days) between successive sampling intervals:
      t.int = diff(t)

      # Do we want to create a wrap around interval for the last interval to complete a year?
      if(wrapDate){
        t.int = c(t.int, (364-sum(t.int)),NA)
        date1 = c(date1, as.Date(date1[length(date1)])+t.int[length(t.int)-1])
      } else{
        date1 = date1[-length(date1)]
      }
      returnDf = list()
      returnDf[[dateCol]] <- as.Date(date1, origin = as.Date("1970-01-01"))
      returnDf[["julianDate"]] <- julian(date1, origin = as.Date("1970-01-01"))
      returnDf[["int_days"]] <- t.int
      returnDf <- as.data.frame(returnDf)
      # check if any of the columns are coercible to date format
    } else if(!any(unlist(lapply(df, function(x) dateCoercible(x))))){
      stop("`dateCol` is not a date or coercible.")
    } else{
      # extract the date column, unlist to vector and maintain class
      dateCol = names(df[unlist(lapply(df, function(x) inherits(x, c('Date', 'POSIXt'))))])
      date1 = unname(do.call("c",unique(df[unlist(lapply(df, function(x) dateCoercible(x)))])))
      # convert to julian dates. Origin is default "1970-01-01"
      julian1 = julian(date1, origin = as.Date("1970-01-01"))

      t = sort(union(julian1, NULL))

      #Create a vector of durations (days) between successive sampling intervals:
      t.int = diff(t)

      # Do we want to create a wrap around interval for the last interval to complete a year?
      if(wrapDate){
        t.int = c(t.int, (364-sum(t.int)),NA)
        date1 = c(date1, as.Date(date1[length(date1)])+t.int[length(t.int)-1])
      }
      returnDf = list()
      returnDf[[dateCol]] <- as.Date(date1, origin = as.Date("1970-01-01"))
      returnDf[["julianDate"]] <- julian(date1, origin = as.Date("1970-01-01"))
      returnDf[["int_days"]] <- t.int
      returnDf <- as.data.frame(returnDf)
    }
  } else{
    if('numeric' %in% class(df[[dateCol]])){
      t = sort(union(df[[dateCol]], NULL))

      #Create a vector of durations (days) between successive sampling intervals:
      t.int = diff(t)

      # Do we want to create a wrap around interval for the last interval to complete a year?
      if(wrapDate){
        t.int = c(t.int, (364-sum(t.int)),NA)
        t = c(t, t[length(t)]+t.int[length(t.int)-1])
      } else{
        t = t[-length(t)]
      }
      returnDf = list()
      returnDf[[dateCol]] <- t
      # returnDf[["julianDate"]] <- julian(date1, origin = as.Date("1970-01-01"))
      returnDf[["int_days"]] <- t.int
      returnDf <- as.data.frame(returnDf)
    } else{
      # extract the date column, unlist to vector and maintain class
      dateCol = names(df[unlist(lapply(df, function(x) inherits(x, c('Date', 'POSIXt'))))])
      date1 = unname(do.call("c",unique(df[unlist(lapply(df, function(x) dateCoercible(x)))])))
      # convert to julian dates. Origin is default "1970-01-01"
      julian1 = julian(date1, origin = as.Date("1970-01-01"))

      t = sort(union(julian1, NULL))

      #Create a vector of durations (days) between successive sampling intervals:
      t.int = diff(t)

      # Do we want to create a wrap around interval for the last interval to complete a year?
      if(wrapDate){
        t.int = c(t.int, (364-sum(t.int)),NA)
        date1 = c(date1, as.Date(date1[length(date1)])+t.int[length(t.int)-1])
      }
      returnDf = list()
      returnDf[[dateCol]] <- as.Date(date1, origin = as.Date("1970-01-01"))
      returnDf[["julianDate"]] <- julian(date1, origin = as.Date("1970-01-01"))
      returnDf[["int_days"]] <- t.int
      returnDf <- as.data.frame(returnDf)
    }
  }
  if(is.null(envData)){
    return(returnDf)
  } else{

    if(wrapDate){
      envCols <-names(envData)[names(envData) %ni% c(dateCol,"julianDate","int_days")]
      # does the environmental data have the same length as the wrapped data
      if(nrow(envData) == length(returnDf[[dateCol]])){
        # if same length, is the last row just NA?
        envNAs <- sapply(envCols, FUN = function(x){
          is.na(envData[nrow(envData), x])
        })
        # if yes, wrap and average the first and last sampling date values.
        if(sum(envNAs) > 0){
          for(i in 1:length(envCols)){
            input <- (envData[evnCols[i], 1] + envData[envCols[i],(nrow(envData)-1)])/2
            envData[envCols[i], nrow(envData)] <- input
          }
          returnDf <- merge(returnDf, envData, by = dateCol, all.x = TRUE)
        }
      } else if(nrow(envData) == (length(returnDf[[dateCol]])-1)){
        envWrap <- c()
        envWrap[[dateCol]] <- returnDf[nrow(returnDf), dateCol]
        for(i in 1:length(envCols)){
          input <- (envData[1,envCols[i]] + envData[(nrow(envData)), envCols[i]])/2
          envWrap[[envCols[i]]] <- input
        }
        returnDf <- merge(returnDf, as.data.frame(envData), by = dateCol, all.x = TRUE)
      }
     }
  }
  ## tests ##
  ### are all dates present in the envData object?
  # if(!all(as.character(unlist(unique(returnDf[[dateCol]])))) %in% as.character(unlist(unique(envData[[dateCol]])))){
  #   stop("Error: all sampling dates are not present in envData object")
  # }
  ## end tests ##
  # returnDf <- merge(returnDf, envData, by = eval(dateCol))
  return(as.data.frame(returnDf))
}
