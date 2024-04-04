#functions to convert body length to mass, 
# and gather and replicate data from wide format of orginal full data set.
gather_body_data <- function(DATA, taxon_info = taxa_info,...){
  taxon_info_i = taxon_info %>%
    filter(SITE == unique(DATA$SITE) & TAXON == unique(DATA$TAXON))
  if(DATA$TAXON == "Sperchon glandulosus"){
    DATA_long <- DATA %>%
      gather(key = 'bodysize', value = 'abundance', 6:dim(DATA)[2]) %>%
      filter(abundance != 0) %>%
      mutate_at(vars(SITE:TAXON), as.character) %>%
      mutate(MASS = (1-(taxon_info_i$LM.p.ash/100))*exp((exp(taxon_info_i$LM.a)+taxon_info_i$LM.b*log(as.numeric(bodysize)))))
    DATA_long <- DATA_long[rep(seq_len(dim(DATA_long)[1]), round(DATA_long$abundance,0)),
                           c('SITE','DATE','HABITAT','TAXON', 'MASS')]
  }
  
  DATA_long <- DATA %>%
    gather(key = 'bodysize', value = 'abundance', 6:dim(DATA)[2]) %>%
    filter(abundance != 0) %>%
    mutate_at(vars(SITE:TAXON), as.character) %>%
    mutate(MASS = (1-(taxon_info_i$LM.p.ash/100))*taxon_info_i$LM.a*(as.numeric(bodysize)^taxon_info_i$LM.b))

  DATA_long <- DATA_long[rep(seq_len(dim(DATA_long)[1]), round(DATA_long$abundance,0)),
                         c('SITE','DATE','HABITAT','TAXON', 'MASS')]
}