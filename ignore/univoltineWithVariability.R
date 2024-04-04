library(tidyverse)
library(stagePop)

# set solver options
solver.options=list(DDEsolver='PBS',tol=1e-4,hbsize=1e5,dt=0.01)


#--------------------- univoltine with variability ----------------------
num.strains = 6
univoltineFunctions <- list(
  reproFunc=function(x,time,species,strain){
    set.seed(42)
    A0 = rnorm(num.strains,75,4); q=rlnorm(num.strains, log(5), 0.1)
    reprod=q[strain]*x$univoltine['adults',strain]*exp(-x$univoltine['adults',strain]/A0[strain])
    return(reprod)
  },
  deathFunc=function(stage,x,time,species,strain){
    #per capita death rate (/d)
    # a1 = (x$univoltine['1mm',1]/1000)
    a=c(0.0068,0.005,0.004,0.003,0.0025,0.002,0.0015,0.0001,1)
    return(a[stage])
  },
  durationFunc=function(stage,x,time,species,strain){
    #duration of each stage in days
    a=c(146,10,18,26,36,40,40,48,1)
    return(a[stage])
  },
  immigrationFunc=function(stage,x,time,species,strain){
    v=0
    if (stage==9){if (time>=0 & time<=1){v=101}}
    return(v)
  },
  emigrationFunc=function(stage,x,time,species,strain){return(0)}
)

#-------------------------------RUN MODEL---------------------------------------

modelOutput = popModel(
  numSpecies=1,
  numStages=9,
  numStrains = 6,
  ICs=list(matrix(0,nrow=9,ncol=6)),
  timeVec=seq(0,3650,1),
  timeDependLoss=FALSE,
  timeDependDuration=FALSE,
  rateFunctions=univoltineFunctions,
  solverOptions=solver.options,
  stageNames=list(c('eggs','0,5mm','1mm','1,5mm','2mm','2,5mm','3mm','3,5mm','adults')),
  speciesNames=c('univoltine'),
  sumOverStrains = FALSE,
  plotStrainsFig = FALSE,
  saveFig=FALSE
)

 stablePop = modelOutput %>% data.frame(check.names = FALSE) %>% dplyr::select(-contains("dot")) %>%
   dplyr::filter(between(time, 365,734)) %>%
   pivot_longer(-time, names_to = 'id', values_to = 'n_m2') %>%
   dplyr::mutate(n_m2 = ifelse(n_m2 < 1,0,n_m2),
                 stage = sapply(strsplit(id,"\\."),"[",2),
                 length = tidyr::replace_na(as.numeric(gsub("mm","", gsub(",",".",stage))),0),
                 mass = 0.0025*(length^2.692)*0.958,
                 strain = sapply(strsplit(id, "\\."),"[",3),
                 date = as.Date(time, origin = as.Date("2022-01-01"))) %>%
   dplyr::select(time,date,strain,stage,length,mass,n_m2)

saveRDS(stablePop, "./ignore/varying_univoltine_full.rds")
stablePop_samp = stablePop %>%
  dplyr::filter(time %in% c(365+(30*seq(0:10)),734))
saveRDS(stablePop_samp, "./ignore/varying_univoltine_samp.rds")
