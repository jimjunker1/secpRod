### Subsetting temperatures from global temp file based on sampling dates
temp_subset = function(t, TEMP){
  int_tempC = matrix(NA, ncol = 1, nrow = length(t)+1)
  if (wrap == T){
    j.wrap = t[length(t)] + (365-sum(t.int))
    t = c(t ,j.wrap)
  }
  stream = as.character(SITE)
  for(i in seq_along(t)){
    temp_sub = subset(TEMP, JULIAN >= t[i] & JULIAN <= t[i+1])
    int_tempC[i] = mean(temp_sub[,as.character(stream)], na.rm = T) 
  }
  df= cbind(t, int_tempC)
  out(df)
}

TEMP = temp_subset(t, TEMP)
