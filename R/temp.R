
binList = lapply(20:3, FUN = function(a){
  breaks = seq(min_size, max_size, length.out = a + 1)
  bins = cut(a, breaks = breaks, include.lowest = TRUE)
  any(table(bins) == 0)
})

which

