create_data_lists = function(DATA,envir = parent.frame(),...) {
  #create a list of sites in data frame
  sites_list <- as.list(unique(DATA$SITE))
  #create list of site data subsets
  #sites_data_list <- map2(list(DATA), sites_list, sites_subset)
  #open the future option once download the developer version of furrr
  sites_data_list <- future_map2(list(DATA), sites_list, sites_subset)
  #create lists of taxa in each site
  taxa_lists <- future_map(sites_data_list, function(x) as.list(unique(x$TAXON)))
  LIST = list(sites_list, sites_data_list, taxa_lists)
  names(LIST) = c("sites_list","sites_data_list", "taxa_lists")
  list2env(LIST, envir = parent.frame())
  #lapply(seq_along(LIST), function(x) {
  #  assign(NAMES[x], LIST[[x]], envir = parent.frame())
  #  })
  #pmap(list(LIST, list(NAMES)), function(x,y) assign(y, x, envir = parent.frame()))
  #return(list(sites_list = sites_list, sites_data_list = sites_data_list, taxa_lists = taxa_lists))
}