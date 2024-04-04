#function to subset and average interval temperature data from streams
#require a file of temperature and samples with julian dates
#set the date file
temp_subset = function(bugs_df, TEMP, wrap = T){
  dates = get.dates(bugs_df, site = unique(as.character(factor(bugs_df$SITE))), habitat = "COBBLE")
  int_tempC = matrix(NA, ncol = 1, nrow = nrow(dates)+1)
  if (wrap == T){
    j.wrap = dates[nrow(dates),'JULIAN'] + (365-(dates$JULIAN[nrow(dates)]-dates$JULIAN[1]))
    dates = c(dates$JULIAN ,j.wrap)
  }
  stream = unique(as.character(factor(bugs_df$SITE)))
  for(i in seq_along(dates)){
    temp_sub = subset(TEMP, JULIAN >= dates[i] & JULIAN <= dates[i+1])
    int_tempC[i] = mean(temp_sub[,as.character(stream)], na.rm = T)
  }
  df= cbind(dates, int_tempC)
  df = df[-dim(df)[1],]
  colnames(df) = c('JULIAN', stream)
  return(data.frame(df))
}