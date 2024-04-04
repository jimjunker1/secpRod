mass_positive = function(bootsdata_wide, bootsdata,...){
  for(k in 5:dim(bootsdata_wide)[2]){
    for(l in 1:nrow(bootsdata_wide)){
      if(bootsdata_wide[l,k] < bootsdata_wide[l,(k-1)]){ 
        date.sub = subset(bootsdata, DATE == as.character(names(bootsdata_wide)[k]))
        x=1
        repeat {
          date.samp = date.sub[sample(1:length(date.sub[,1]),1, replace = TRUE),]
          x = x+1
          if(date.samp[,'MASS'] >= bootsdata_wide[l,k-1] | x == 10000){
            bootsdata_wide[l,k] <- date.samp[,'MASS']
            break
          }
        }
      }
    }
  } 
  return(bootsdata_wide) 
}