date_order_lists = function(site_data_list,...){
  date_order = site_data_list %>%
      group_by(DATE) %>%
      dplyr::summarise(mean_mass = median(MASS, na.rm=T), julian = mean(julian)) %>%
      arrange(mean_mass) %>%
      mutate(DATE = reorder(DATE, mean_mass)) %>%
      mutate(day = c(diff(julian),NA), id = 1:n()) %>%
      mutate(day = ifelse(day < 0, 365-abs(day),day))
}