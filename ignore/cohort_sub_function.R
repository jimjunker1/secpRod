######  Function work #####
cohort_subset = function(cohort_data, sample_data = hengill_full, parallel = FALSE){
  library(plyr)
  library(dplyr)
  path = "./Output/stream_cohorts/"
  cohort_data <- cohort_data %>%
    mutate_at(vars(SITE:COHORT), as.factor)
  #sample_data <- sample_data %>%
  #  mutate_at(vars(SITE:TAXON),as.factor)
  if(parallel == FALSE){
  # cohort_data[,"SITE":"COHORT"] = lapply(cohort_data[,SITE:COHORT], factor)#;col_names = names(PROP_TABLE)
  # #sample_data[,"SITE":"COHORT"] = lapply(sample_data[,SITE:COHORT], factor)#;col_names = names(PROP_TABLE)
  # sample_data = sample_data
  sites = droplevels(cohort_data$SITE)
for(i in sites){
    data_site = cohort_data[which(cohort_data$SITE == i),]
    taxon = droplevels(data_site$TAXON)
    
    for(j in taxon){
      data_taxon = data_site[which(data_site$TAXON == j),]
      cohorts = droplevels(data_taxon$COHORT)
      
      #subset of taxon samples
      #sample_sub = subset(sample_data, SITE == as.character(unique(data_taxon$SITE)) & TAXON == as.character(unique(data_taxon$TAXON)))
      #sample_sub = subset(sample_data, SITE == as.character(i) & TAXON == as.character(j))
      sample_sub = sample_data[which(sample_data$SITE == i & sample_data$TAXON == j),]
      #sample_sub = sample_data[which(sample_data$SITE == as.character(i) & sample_data$TAXON == as.character(j)),]
      #sample_sub = sample_data[which(sample_data$SITE == i),]
      
      for(k in cohorts){
        data_cohort = data_taxon[which(data_taxon$COHORT == k),]
        dates = droplevels(data_cohort$DATE)
        
        newdata = data.frame(matrix(ncol = dim(sample_data)[2], nrow = 0))
        colnames(newdata) = names(sample_data)
        for(l in dates){
          date_sub = subset(sample_sub, sample_sub$DATE == l)# sample data of date
          cohort_date = data_cohort[which(data_cohort$DATE == l),]
          cohort_sizes = date_sub[,6:dim(date_sub)[2]]#just the size class columns
          size_cols = date_sub[,1:5]#just the header columns
          size_sub = cohort_sizes[,which(as.numeric(colnames(cohort_sizes)) >= as.numeric(cohort_date$START_SIZE) &
                                           as.numeric(colnames(cohort_sizes)) <= as.numeric(cohort_date$END_SIZE))]
          cohort_full = cbind(size_cols,size_sub)
          newdata = bind_rows(newdata, cohort_full)
         }
        newdata[is.na(newdata)] <- 0
        #print(newdata)
        #print(paste(i,"_",j,"_",k))
        #write.csv(newdata, file = paste(as.character(i),"_",as.character(j),"_",as.character(k),".csv", sep = ""), row.names = F)
        write.table(newdata, file = paste(path, as.character(i),"_",as.character(j),"_",as.character(k),".txt", sep = ""), sep = "\t", quote = F, row.names = F)
      }
    }
  }  
  } else{
    cohort_list = cohort_data %>%
      group_by(SITE, TAXON, COHORT) %>%
      do(data = (.)) %>% pull(data)
    
    cohort_data_list <<- map(cohort_list, ~.x %>%
                             pmap(~data.frame(sample_data %>%
                                                mutate_at(vars(SITE:TAXON),as.factor) %>%
                                                filter(SITE == ..1 & DATE == ..2 & TAXON == ..4) %>%
                                                select(SITE:TAXON, which(as.numeric(colnames(.)) >= ..6 & as.numeric(colnames(.)) <= ..7)),
                                              check.names = F)))
    
    # data_subset = function(cohort_row, sample_data = hengill_full,...){
    #   pmap(cohort_row, function(x){
    #     site_subset = sample_data %>%
    #       filter(SITE == unique(cohort_row$SITE) & TAXON == unique(cohort_row$TAXON) & DATE == droplevels(cohort_row$DATE)) %>%
    #       select_if(function(col) as.numeric(col) >= cohort_row$START_SIZE &
    #                   as.numeric(col) <= cohort_row$END_SIZE)})}
    #  
      #for(i in 1:length(droplevels(cohort_data$DATE))){ 
      #cohort_sizes = date_sub[,6:dim(date_sub)[2]]#just the size class columns
      #size_cols = date_sub[,1:5]#just the header columns
      #cohort_full = sample_data %>%
       # filter(DATE == i) %>%
        #select_if(function(col) as.numeric(col) >= data_subset$START_SIZE &
        #            as.numeric(col) <= data_subset$END_SIZE)
      
      #size_sub = cohort_sizes[,which(as.numeric(colnames(cohort_sizes)) >= as.numeric(cohort_date$START_SIZE) &
      #                                 as.numeric(colnames(cohort_sizes)) <= as.numeric(cohort_date$END_SIZE))]
      #cohort_full = cbind(size_cols,size_sub)
      #newdata = bind_rows(newdata, cohort_full)
     
      #cohort_data_list = map(cohort_list, data_subset) 
  
      #return(cohort_data_list = cohort_data_list)
    }
    
    }
#  }
#}
