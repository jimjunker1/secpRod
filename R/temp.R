x = rgamma(1e3, 1,2);hist(x)
min_size = min(x)
max_size = max(x)
binList = sapply(20:3, FUN = function(a){
  breaks = seq(min_size, max_size, length.out = a + 1)
  bins = cut(x, breaks = breaks, include.lowest = TRUE)
  zeros = any(unlist(table(bins)) == 0)
  return(zeros)
})

which

a = 20
!binList
