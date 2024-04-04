#boot strap function

parallel_cohort_boot = function(DATA, tax_info = taxa_info,nboot = NULL, parallel = TRUE,...){
  library(tidyverse)
  library(furrr)
  source("./Analysis/functions/cohort_boot_functions/gather_data_function.R")
  source("./Analysis/functions/cohort_boot_functions/create_data_lists_function.R")
  source("./Analysis/functions/cohort_boot_functions/sites_subset_function.R")
  source("./Analysis/functions/cohort_boot_functions/convert_to_julian_function.R")
  source("./Analysis/functions/cohort_boot_functions/taxa_list_split_function.R")
  source("./Analysis/functions/cohort_boot_functions/date_order_lists_function.R")
  source("./Analysis/functions/cohort_boot_functions/date_reorder_function.R")
  source("./Analysis/functions/cohort_boot_functions/run_list_boots_function.R")
  source("./Analysis/functions/cohort_boot_functions/mass_positive_function.R")
  source("./Analysis/functions/cohort_boot_functions/calculate_growth_function.R")
  source("./Analysis/functions/cohort_boot_functions/join_days_function.R")

  #sets the number of "individuals" you want to sample
  if(is.null(nboot)){
    nboot = 5
  } else{ nboot = nboot }
  if(parallel == TRUE) {
    set.seed(100) # sets randomization seed for bootstrapping, turn on if you want the same results every time you run the code 
    #DATA <- DATA %>% select(-JULIAN)
    DATA = gather_data(DATA, tax_info)
    #plan(multiprocess)
    create_data_lists(DATA)
    #debugonce(convert_to_julian)
    sites_data_list = map(sites_data_list, convert_to_julian)
    #site_taxa_data_lists = future_pmap(list(sites_data_list,taxa_lists), taxa_list_split)
    cohort_date_lists = map(sites_data_list, date_order_lists)
    site_taxa_data_lists = date_reorder(sites_data_list,cohort_date_lists)
    source("./Analysis/functions/cohort_boot_functions/parallel_bootstrap_function.R")#function selects data 
    #debugonce(boots)
    #debugonce(date_order_lists);
    #list of lists of bootstapped body sizes for each site-taxa/cohort
    #bootsdata = future_map2(site_taxa_data_lists,nboot, run_list_boots)
    # browser();debugonce(run_list_boots)
    bootsdata = run_list_boots(site_taxa_data_lists, nboot)
    #now work across two columns at a time to estimate growth
    #1:names(bootsdata[,-1:2]))
    #browser();
    bootsdata = join_days(bootsdata, cohort_date_lists)
    # debugonce(calculate_growth)
    bootsdata = calculate_growth(bootsdata, cohort_date_lists)
    bootsdata = map2(bootsdata, cohort_date_lists, ~.x %>% left_join(.y %>% select(start_date = DATE, mean_mass)))
    
    TAXON = as.character(unique(droplevels(as.factor(flatten(bootsdata)$TAXON))))
    SITE = as.character(unique(droplevels(as.factor(flatten(bootsdata)$SITE))))
    #names(bootsdata) = sites_list
    #future_map2(sites_list, taxa_lists, function(x,y) {
    #  future_map2_chr(x, y, function(site, taxon) paste(site,"_",taxon,sep=""))})
    #setNames(lapply(bootsdata, setNames, taxa_lists), sites_list)
    write.csv(bootsdata, file = paste("./Output/",TAXON,"_site-",SITE,"_IGR.csv",sep = ""), row.names = F) 
    return(bootsdata)
 } else{
    
  taxa = list(unique(levels(DATA$species)))#set taxa levels
  tax_data = map(taxa, taxa_subset)
  
  for(i in taxa){
  tax_data = DATA[which(DATA$species == as.character(i)),]#grab single taxa
  sites = unique(levels(as.factor(tax_data$Site)))#set site levels for a single taxa
for(j in sites){
  site_data = tax_data[which(as.factor(tax_data$Site) == j),]#grab all data for a single site

#convert to dates and calculate julian day
site_data$Pd = as.POSIXct(site_data$date, format = "%d-%m-%y")
site_data$date = as.factor(site_data$Pd)
site_data$julian = as.numeric(format(site_data$Pd, "%j"));site_data$Pd = NULL

#reorder dates based start of cohort
site_data %>%
  group_by(date) %>%
  summarize(mean_mass = mean(mass, na.rm=T), julian = unique(julian, na.rm = T)) %>%
  arrange(mean_mass) %>%
  mutate(date = reorder(date, mean_mass)) %>%
  mutate(day = c(diff(julian),NA), id = 1:n()) -> date_df

#fix the weird wrap around dates
wrap_loc = which(date_df$day < 0);wrap_fix = 365-abs(date_df[which(date_df$day < 0),'day'])
date_df[wrap_loc, 'day'] = wrap_fix

#reorder the site data based on 'cohort date'
site_data$date = factor(site_data$date, levels = levels(date_df$date))

source("./bootstrap_function.R")#function selects data 
#debugonce(boots)
bootsdata = boots(site_data, nboot = nboot)

bootsdata %>%
  select(species, Site, date, mass) %>%
  left_join(date_df[,c('date','id')]) -> bootsdata_id

bootsdata %>%
  select(species, Site, date, mass) %>%
  mutate(id = rep(1:nboot, nrow(date_df))) %>%
  spread(date,mass) %>%
  select(-id)-> bootsdata_wide

#now work across two columns at a time to estimate growth
#1:names(bootsdata[,-1:2]))
growth_rate <- function(x,y,z) {(log(y/x))/z}

igr_df = c()
for(w in 1:(max(bootsdata_id$id)-1)){
  igr_vec = mapply(growth_rate, bootsdata_id[which(bootsdata_id$id == w),'mass'], bootsdata_id[which(bootsdata_id$id == w+1),'mass'],
                   date_df[w,'day'])
  igr_df = cbind(igr_df,igr_vec)
}

colnames(igr_df) = unique(levels(date_df$date))[1:(nrow(date_df)-1)]
igr_df = data.frame(species = as.character(i), site = j, igr_df, check.names = F)
igr_df %>%
  gather(start_date, IGR, 3:(dim(igr_df)[2]), factor_key = T) -> igr_long

igr_fix = which(igr_long$IGR <= 0)
igr_long[igr_fix, 'IGR'] = 0.001

write.csv(igr_long, file = paste("./output/",i,"_site-",j,"_IGR.csv",sep = ""), row.names = F) 
  }  
 }  
}
} 