convert_to_julian = function(data_list,...) {
  #posix_vec = as.POSIXct(data_list$DATE, format = "%m/%d/%y")
  posix_vec = as.POSIXct(data_list$DATE)
  data_list %>% mutate(DATE = as.factor(posix_vec)) %>%
    mutate(julian = as.numeric(format(posix_vec, "%j"))) %>%
    mutate_at(vars(DATE), as.character)
}