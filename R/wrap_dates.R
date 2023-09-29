#'
#'
#'
wrap_dates <- function(df = NULL, dateCol = NULL, wrapDate = TRUE,...) {
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
      return(as.data.frame(returnDf))
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
      return(as.data.frame(returnDf))
    }
  } else{

  }
}
