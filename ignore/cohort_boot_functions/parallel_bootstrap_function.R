# Bootstrap function
parallel_boots = function(allData, nboot,...){
  suppressMessages(library(parallel))
  suppressMessages(library(doParallel))
  suppressMessages(library(foreach))
  # for each date resample data
  boot.data = c()
  for(j in unique(levels(allData$DATE))){ # for each date
    date.sub <- subset(allData, DATE == j)
    #resample rows #s of date.sub
    date.samp <- sample( 1:length(date.sub[,1]), size=nboot, replace=TRUE )
    #compile resampled data
    date.samp <- date.sub[date.samp,]
    boot.data = rbind(boot.data,date.samp)
  }
  boot.data %>%
    select(TAXON, SITE, DATE, start_size, end_size, MASS) %>%
    mutate(id = rep(1:nboot, length(unique(levels(allData$DATE))))) %>%
    spread(DATE,MASS) %>%
    select(-id) -> boot.data_wide

  # debugonce(mass_positive)
    boot.data_wide = mass_positive(boot.data_wide, boot.data)
    boot.data_wide %>%
      group_by(TAXON, SITE) %>%
      gather(key = DATE, value = MASS, 5:dim(boot.data_wide)[2])-> boot.data
  
  ### set up the doParallel foreach loop ##
  # make sure mass change is positive
  #foreach(date_column in 4:dim(boot.data_wide)[2]) %dopar% {
  #  foreach(row_i in 1:nrow(boot.data_wid)) %dopar% {
  #    if(boot.data_Wide[row_i,date_column] < boot.data_wide[row_i,(date_column + 1)]){
  #      date.sub = subset(allDATE, DATE == as.character(names(boot.data_wide)[date_column]))
  # tic()
  # for(k in 4:dim(boot.data_wide)[2]){
  #   for(l in 1:nrow(boot.data_wide)){
  #     if(boot.data_wide[l,k] < boot.data_wide[l,(k-1)]){ 
  #       date.sub = subset(allData, DATE == as.character(names(boot.data_wide)[k]))
  #       x=1
  #         repeat {
  #       date.samp = date.sub[sample(1:length(date.sub[,1]),1, replace = TRUE),]
  #       x = x+1
  #       if(date.samp[,'MASS'] >= boot.data_wide[l,k-1] | x == 300){
  #         boot.data_wide[l,k] <- date.samp[,'MASS']
  #         break
  #         }}
  #     }
  #   }
  #   
  # };toc()
  #clean up and output data
  rownames(boot.data) = NULL
  #list2env(boot.data, envir = parent.frame())
  return (boot.data)
}