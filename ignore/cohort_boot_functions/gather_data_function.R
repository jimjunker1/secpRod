#functions to convert body length to mass, 
# and gather and replicate data from wide format of orginal full data set.
gather_data <- function(DATA, taxon_info = taxa_info,...){
  taxon_info_i = taxon_info %>%
    filter(SITE == unique(DATA$SITE) & TAXON == unique(DATA$TAXON))
  
  DATA_long <- DATA %>%
    gather(key = 'bodysize', value = 'abundance', 9:dim(DATA)[2]) %>%
    filter(abundance != 0) %>%
    mutate_at(vars(start_size:TAXON), as.character) %>%
    mutate(MASS = (1-(taxon_info_i$LM.p.ash/100))*taxon_info_i$LM.a*(as.numeric(bodysize)^taxon_info_i$LM.b))

  DATA_long <- DATA_long[rep(seq_len(dim(DATA_long)[1]), round(DATA_long$abundance,0)),
                         c('SITE','DATE','HABITAT','TAXON','start_size','end_size','MASS')]
}