join_days = function(site_taxa_data_list, cohort_date_list,...) {
  future_map2(site_taxa_data_list, cohort_date_list, function(x,y) {
    x %>% left_join(y %>% select(DATE, day, id))
  })
}