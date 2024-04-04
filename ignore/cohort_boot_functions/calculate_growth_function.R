calculate_growth = function(bootdata, cohort_date_list,...){
  growth_rate <- function(x,y,z) {(log(y/x))/z}
  future_map2(bootdata, cohort_date_list, function(x,y){
  # map2(bootdata, cohort_date_list, function(x,y){
    igr_df = c()
    for(date_id in 1:(max(x$id)-1)){
      igr_vec = mapply(growth_rate, x[which(x$id == date_id),'MASS'], x[which(x$id == date_id+1),'MASS'],
                       y[date_id,'day'])
      igr_df = cbind(igr_df,igr_vec)}
  colnames(igr_df) = unique(levels(y$DATE))[1:(nrow(y)-1)]
  igr_df = data.frame(TAXON = as.character(unique(x$TAXON)), SITE = as.character(unique(x$SITE)),
                      start_size = unique(x$start_size), end_size = unique(x$end_size), igr_df, check.names = F)
  #View(igr_df)
  igr_df %>%
    gather(start_date, IGR, 5:(dim(igr_df)[2]), factor_key = T) -> igr_long
  igr_fix = which(igr_long$IGR <= 0)
  igr_long[igr_fix, 'IGR'] = 0.001
  # TAXON = as.character(unique(droplevels(as.factor(x$TAXON))))
  # SITE = as.character(unique(droplevels(as.factor(x$SITE))))
  #igr_long = future_map2(igr_long, cohort_date_lists, ~.x %>% left_join(.y %>% select(start_date = DATE, mean_mass)))
  #write.csv(igr_long, file = paste("./Output/",TAXON,"_site-",SITE,"_IGR.csv",sep = ""), row.names = F) 
  return(igr_long)})
}