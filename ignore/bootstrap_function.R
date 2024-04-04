# Bootstrap function
boots = function(allData){
  # for each sampling day resample data
  boot.data = 0.
  for (j in 1:length(levels(as.factor(allData$DATE)))){ # for each date
    
    date.sub <- subset(allData, DATE == levels(as.factor(allData$DATE))[j])
    #resample rows #s of date.sub
    date.samp <- sample( 1:length(date.sub[,1]), size=nReplicates, replace=TRUE )
    #compile resampled data
    date.samp <- date.sub[date.samp,]
    if (j==1) {
      boot.data = date.samp
    } else {
      boot.data = rbind(boot.data,date.samp)
    }
  }
  #clean up and output data
  rownames(boot.data) = NULL
  return (boot.data)
}