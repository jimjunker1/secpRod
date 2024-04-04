# population modeling of different lifehistories
library(tidyverse)
library(stagePop)
file.edit(system.file("DemoFiles/MultipleStrains.R", package = 'stagePop'))

### estimate population structure of univoltine population with fixed duration egg stage
proj.mat = function(n0, matproj, tmax){
  res.mat<-matrix(NA,nrow=tmax+1,ncol=length(n0))
  res.mat[1,]<-n0
  for(i in 2:(tmax+1))
  {
    nA = res.mat[(i-1),9]
    res.mat[i,]=matproj %*% res.mat[(i-1),]
  }
  return(res.mat)
}


# create life table matrix
nA = 10
A = matrix(rbind(c(0,0,0,0,0,0,0,0,(5*nA*exp(-nA/20))),
                 c(0.8,0,0,0,0,0,0,0,0),
                 c(0,0.8,0,0,0,0,0,0,0),
                 c(0,0,0.8,0,0,0,0,0,0),
                 c(0,0,0,0.2,0,0,0,0,0),
                 c(0,0,0,0,0.6,0,0,0,0),
                 c(0,0,0,0,0,0.7,0,0,0),
                 c(0,0,0,0,0,0,0.9,0,0),
                 c(0,0,0,0,0,0,0,0.9,0)), nrow = 9, ncol = 9)

N0 = matrix(c(50,40,32,25,20,3,2,2,2), ncol = 1)

nEst = proj.mat(n0 = N0, matproj = A, tmax = 2555)

nPop<-apply(nEst,1, sum)

propEst<-nEst/nPop
matplot( propEst, type="l")

popEst = nEst[366:731,5:8]
matplot(popEst,type = 'l')
plot(nPop)





## using stagePop to estimate population with egg stage

solver.options=list(DDEsolver='PBS',tol=1e-4,hbsize=1e5,dt=0.01)

#--------------------------DEFINE RATE FUNCTIONS---------------------------------------
univoltineFunctions <- list(
  reproFunc=function(x,time,species,strain){
    A0 = 75; q=5
    reprod=q*x$univoltine['adults',1]*exp(-x$univoltine['adults',1]/A0)
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
  ICs=list(matrix(0,nrow=9,ncol=1)),
  timeVec=seq(0,3650,1),
  timeDependLoss=FALSE,
  timeDependDuration=FALSE,
  rateFunctions=univoltineFunctions,
  solverOptions=solver.options,
  stageNames=list(c('eggs','0.5mm','1mm','1.5mm','2mm','2.5mm','3mm','3.5mm','adults')),
  speciesNames=c('univoltine'),
  saveFig=FALSE
)

stablePop = modelOutput[,1:9]
stablePop[stablePop<1] <- 0
univoltine = stablePop[365:735,]
saveRDS(univoltine, "./ignore/univoltine_full.rds")
saveRDS(univoltine_samp, "./ignore/univoltine_samp.rds")

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
  stageNames=list(c('eggs','0.5mm','1mm','1.5mm','2mm','2.5mm','3mm','3.5mm','adults')),
  speciesNames=c('univoltine'),
  sumOverStrains = FALSE,
  plotStrainsFig = FALSE,
  saveFig=FALSE
)

stablePop = modelOutput[,1:9]
stablePop[stablePop<1] <- 0
univoltine = stablePop[365:735,]
saveRDS(univoltine, "./ignore/varying_univoltine_full.rds")
saveRDS(univoltine_samp, "./ignore/varying_univoltine_samp.rds")

## semivoltine
#--------------------------DEFINE RATE FUNCTIONS---------------------------------------
multivoltineFunctions <- list(
  reproFunc=function(x,time,species,strain){
    A0 = 75; q=10
    reprod=q*x$multivoltine['adults',1]*exp(-x$multivoltine['adults',1]/A0)
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
    a=c(5,1,3,4,5,5,5,5,5)
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
  ICs=list(matrix(0,nrow=9,ncol=1)),
  timeVec=seq(0,3650,1),
  timeDependLoss=FALSE,
  timeDependDuration=FALSE,
  rateFunctions=multivoltineFunctions,
  solverOptions=solver.options,
  stageNames=list(c('eggs','0.5mm','1mm','1.5mm','2mm','2.5mm','3mm','3.5mm','adults')),
  speciesNames=c('multivoltine'),
  saveFig=FALSE
)

stablePop = modelOutput[,1:9]
stablePop[stablePop<1] <- 0
multivoltine = stablePop[365:735,]
saveRDS(multivoltine, "./ignore/multivoltine_full.rds")
multivoltine_samp = multivoltine %>% data.frame %>% dplyr::mutate(date = as.Date(time, origin = as.Date("2022-01-01"))) %>% dplyr::select(date, doy = 'time', everything()) %>%
  .[c(2,32,62,92,122,152,182,212,242,272,302,332),]
saveRDS(multivoltine_samp, "./ignore/multivoltine_samp.rds")
