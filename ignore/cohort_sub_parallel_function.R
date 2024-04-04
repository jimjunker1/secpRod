######  Function work #####
cohort_subset_parallel = function(cohort_data, sample_data, parallel = FALSE){
if(parallel == TRUE){
  plan(multiprocess)
  cohort_list = cohort_data %>%
    group_by(SITE, TAXON, COHORT) %>%
    do(data = (.)) %>% pull(data)
  
  cohort_data_list = future_map(cohort_list, ~.x %>%
                             future_pmap(~data.frame(sample_data %>%
                                                mutate_at(vars(SITE:TAXON), as.factor) %>%
                                                filter(SITE == ..1& DATE == ..2& TAXON == ..4) %>%
                                                select(SITE:TAXON, which(as.numeric(colnames(.)) >= ..6 & as.numeric(colnames(.)) <= ..7)),
                                              check.names = F)))
  return(cohort_data_list)
} else {
  cohort_list = cohort_data %>%
    group_by(SITE, TAXON, COHORT) %>%
    do(data = (.)) %>% pull(data)
  
  cohort_data_list = map(cohort_list, ~.x %>%
                                  pmap(~data.frame(sample_data %>%
                                                            #mutate_at(vars(SITE:TAXON,-JULIAN), as.factor) %>%
                                                            filter(SITE == ..1& DATE == ..2& TAXON == ..4) %>%
                                                            select(SITE:TAXON, which(as.numeric(colnames(.)) >= ..6 & as.numeric(colnames(.)) <= ..7)),
                                                          check.names = F)))
}
}