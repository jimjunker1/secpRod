#2oP Bootstrap, Author: Cam MacKenzie; Modified by: Jim Junker
#cmackenzie@atwaterresources.com
#james.junker1@gmail.com
#last modified: April 2018
#Coded from 2oP macro 
#========================================================
# Input data 
#==============================================================================
# Data (.txt) MUST be in the following format, in particular the date (i.e. %m/%d/%y), 
# but number of columns before size class data should also be the same:
# DATE    SAMPLE  SITE  HABITAT  TAXON    1 2 3 4 5 6.....
# 6/18/12	1       Hver  COBBLE   Baetidae	0	0	0	0	0	0	0
# 6/18/12	2       Hver  COBBLE   Baetidae	44	0	0	0	0	0
# 6/18/12	3	      Hver  COBBLE   Baetidae	0	0	0	0	0	0
# 6/18/12	4	      Hver  COBBLE   Baetidae	0	0	0	0	0	0
# 6/18/12	5	      Hver  COBBLE   Baetidae	44	0	0	0	0	0
# 6/18/12	6	      Hver  COBBLE   Baetidae	0	0	0	0	0	0
# 6/27/12 1	      Hver  COBBLE   Baetidae	77	444	0	0	0	0
# 6/27/12	2	      Hver  COBBLE   Baetidae	266	44	0	0	0	0
# 6/27/12	5	      Hver  COBBLE   Baetidae	79	88	0	0	0	0
#...

# Output format
#===============================================================================
# Data (.txt) will be written with the following columns
# empty   boot.mean.1 2.50% 97.50%
# N <- mean abundance of all samples-dates
# B <- mean biomass of all samples-dates
# P <- mean production estimate across all sample-dates
# igr <- instantaneous growth rate over all samples-dates
# P_B <- P/B ratio estimate for all samples-dates
# Pint.1 <- Production estimate over first interval
# Pint.2 <- Production estimates over second interval
# .
# .
# Pint.n <- Production estimates over the nth (last) interval
# Bint.1 <- Biomass estimates over the first interval
# .
# .
# Bint.n <- Biomass estimates over the nth (last) interval
# Nint.1...Nint.n <- Abundance estimates over 1st...nth interval
# size.1...size.n <- mean mass estimate over 1st...nth interval
# Pd1...Pdn <- Daily production flux over 1st...nth interval
# Pd.B1...Pd.Bn <- instaneous growth rate for 1st...nth interval (d^-1)
#######
# igr.boot <- function(DATA, TAXA,...)

#Input and Output file names
#Crop specific dates for use in IGR model dates
#Be sure individual mass values are appropriate for species
#  size classes used
#Format Date
#Set appropriate # of Production Intervals

#graphics.off()
#rm(list=ls())
set.seed(100) # sets randomization seed for bootstrapping, turn on if you want the same results every time you run the code 

IGR_file_boot <- function(path, file.names,tax.info, nBoot){
  library(data.table)
  for(a in file.names){
  allData <- read.table(file = paste(path,a, sep = ""),header=TRUE, sep = "\t",check.names = F,stringsAsFactors=FALSE)
  file.name.out = paste(path,"/IGR_boot_output/",gsub(x = a, '.txt', replacement = ""),"_prod.csv", sep = "") 
  
  #remove any sneaky rows with all zeros
  allData = allData[rowSums(allData[,6:dim(allData)[2]])>0,]
  # Define samples to crop (due to outliers) for igr calculation  
  # enter 0 (don't crop) or crop samples as follows c(1,9) or 9;   
  # removes first and ninth sample or ninth sample only, respectively
  # for igr calculation 
  crop.dates = 0
  # Define number of Monte Carlo simulations to run
  nBoot= nBoot
  print(a)
  # Define number replicates for the bootstrap (usually 6). 
  # This variable alters the number of replicates the bootstrap
  # produces per sampling date
  nReplicates<<-5 #adjust this eventually to set automatically based on actual samples
  
#Create a dataframe of sizes:
  #print(unique(allData$TAXON))
  stream = unique(allData$SITE)
  cohort = substring(a, nchar(a) - 4, nchar(a) -4)
  taxon = tax_info[which(tax_info$TAXON == unique(allData$TAXON)),1:5]
  mm <- as.numeric(substr(names(allData)[6:(dim(allData)[2])], 1, 6))
  AFDM <- (1-(taxon$LM.p.ash/100))*taxon$LM.a*(mm^taxon$LM.b)
  sizes <- data.frame(mm, AFDM)
  
  ind.mass = sizes$AFDM
  #ind.size = sizes$mm#c(0.00003, 0.00025, 0.00084, 0.00200, 0.00392, 0.00678, 0.01078, 0.01612, 0.02299, 0.03157, 0.04206, 0.05466, 0.06955, 0.08694,	0.10701, 0.12997, 0.15599, 0.18529, 0.21805, 0.25447, 0.29473, 0.33905,	0.38761, 0.44060, 0.49822, 0.56068, 0.62815, 0.70084, 0.77895, 0.86266,	0.95218, 1.04769, 1.14940, 1.25751, 1.37220, 1.49367, 1.62212, 1.75774, 1.90074, 2.05130, 2.20963, 2.37592, 2.55036, 2.73315) #mass for each size class
#source("./Analysis/functions/bootstrap_function.R")  
source("./Analysis/functions/cohort_boot_functions/parallel_bootstrap_function.R")
  #browser()
  # Monte Carlo simulations
  #==============================================================================
  # To test calculations compared to excel macro use:
  # bootsData = read.csv("2oP_boots_data.csv", header = TRUE) 
  #MCsims = data.frame(matrix(NA ,nrow = nBoot, ncol = 14))
  #MCsum = data.frame(nrow = nBoot, ncol = 8)
  #MCsims = matrix(ncol = 14)
  MCsims = c()
  #MCsims = matrix(nrow = nBoot * (length(unique(bootsData$DATE))-1), ncol = 14)
  MCsum = matrix(nrow = nBoot, ncol = 8)
  for (n in 1:nBoot){
    
    bootsData = boots(allData)
    
    size.classes = bootsData[,6:length(bootsData)] #extract size class data
    sum.classes = colSums(size.classes) # sum columns
    
    # determine maximum size class, last column with data
    no.classes = max(which(sum.classes != 0))
    
    #for (j in length(sum.classes):1){
    #  if (sum.classes[j]>0){
    #   no.classes = j
    #    break
    #  }
    #}
    
    #trim data - remove extra columns with no data
    bootsData = bootsData[,1:(no.classes+5)]
    
    #format dates
    bootsData$DATE = as.Date(as.character(bootsData$DATE), format = "%m/%d/%y")
    
    # Calculate Pararmeters N,B,P,igr
    #==============================================================================
    #no.classes = number size classes = k
    dates = levels(as.factor(bootsData$DATE)) #  = i
    
    ## 1) Calculate Abundance (N) table
    N.tab = matrix( ncol = no.classes+3, nrow = length(dates))
    
    #for (i in 1:length(dates)){
    #  temp.dat = 0.
    #  date.sub = subset(bootsData, DATE == dates[i])
    #  temp.dat[1] = dates[i]
    #  for (k in 1:no.classes){
    #    #calculate mean of bootstrapped size class for date i, rounded as in macro
    #    temp.dat[k+1] = as.numeric(sum(date.sub[k+5]))/length(date.sub[,1]) #values are rounded as in macros
    #  }
    #  N.tab = rbind(N.tab,temp.dat)
    #}
    
    #N.tab = N.tab[-1,]
    #rownames(N.tab) = NULL
    #colnames(N.tab) = c("Date", 1:no.classes)
    #N.tab[,1] = as.Date(N.tab[,1]) 
    
    for (i in 1:length(dates)){
      date.sub = subset(bootsData, DATE == dates[i])
      N.tab[i,1] = dates[i]
      size.classes = date.sub[,6:length(date.sub)] #extract size class data
      sum.classes = colSums(size.classes) # sum columns
      
      if(is.na(max(which(sum.classes != 0))) | is.infinite(max(which(sum.classes != 0)))){
        next()
      } else(j = max(which(sum.classes != 0)))
      
      END_SIZE = as.numeric(colnames(size.classes[as.numeric(j)]))
      NonzeroIndex = as.numeric(min(which(sum.classes > 0)))
      if(is.null(NonzeroIndex) | is.na(NonzeroIndex) | is.infinite(NonzeroIndex)){
        START_SIZE = 0.25
      } else(START_SIZE = as.numeric(colnames(size.classes[as.numeric(NonzeroIndex)])))
      
      for (k in 1:no.classes){
        #calculate mean of bootstrapped size class for date i, rounded as in macro
        #N.tab[i, k+1] = as.numeric(sum(date.sub[i,k+5]))/length(date.sub[i,1]) #values are rounded as in macros--old
        N.tab[i,k+1] = as.numeric(sum(date.sub[,k+5]))/nReplicates #values are rounded as in macros
      }
      
      N.tab[i, no.classes+2] = START_SIZE
      N.tab[i, no.classes+3] = END_SIZE
    }
    rownames(N.tab) = NULL
    colnames(N.tab) = c("Date", 1:no.classes, "START_SIZE", "END_SIZE")
    N.tab[,1] = as.Date(N.tab[,1]) 
    ## 2) Calculate Biomass (B) table
   #need to quite dealing with this vector of varying types crap
    Biomass.tab = matrix(ncol = no.classes+6, nrow = length(dates))
    for (i in 1:length(dates)){
      Biomass.tab[i,1] = dates[i]
      Biomass.tab[i,2] = sum(as.numeric(N.tab[i,1:no.classes+1]))

      for (k in 1:no.classes){
        Biomass.tab[i,k+2] = as.numeric(N.tab[i,k+1])*ind.mass[k]
      }
      Biomass.tab[i,no.classes+3] = sum(as.numeric(Biomass.tab[i,1:no.classes+2]))/as.numeric(Biomass.tab[i,2]) #calculate avg mass
      Biomass.tab[i,no.classes+4] = sum(as.numeric(Biomass.tab[i,1:no.classes+2])) #calculate B mg/m2
      Biomass.tab[i,no.classes+5] = N.tab[i,no.classes+2]
      Biomass.tab[i,no.classes+6] = N.tab[i,no.classes+3]
      }
    
    rownames(Biomass.tab) = NULL
    colnames(Biomass.tab) = c("Date", "Abundance", 1:no.classes, "avg mass", "Biomass", "START_SIZE", "END_SIZE") 
    
    ## 3) Estimate igr 
    
    #define D.elapse
    dates = as.Date(dates, format = "%Y-%m-%d")
    
    D.elapse = 0.
    for (j in 1:(length(dates)-1)){
      D.elapse[j] = dates[j+1]-dates[1]
    }
    D.elapse = c(0,D.elapse)
    
    
    #log transform avg mass
    #browser()
    ln.avg.mass = log(as.numeric(as.character(Biomass.tab[,"avg mass"])))
    #print(ln.avg.mass)
    #Crop dates for appropriate IGR model (specified by crop.dates above)
    if (crop.dates[1] != 0){
      D.elapse = D.elapse[-crop.dates]
      ln.avg.mass = ln.avg.mass[-crop.dates]
    }
    
    # 
    if(all(is.na(ln.avg.mass)) | all(is.na(D.elapse))){
      next()
    }
    #calculate igr from change in mass
    igr = (ln.avg.mass[2]-ln.avg.mass[1])/D.elapse[2]#;print(igr)
    
###### old code. Don't need to run lm based on this method I am using#####   
     # Run linear model
    #growth.lm = lm(ln.avg.mass ~ D.elapse)
    
    #summary.lm(growth.lm)
  #browser()
    # Diagnostic plots, if wanted
    #plot(growth.lm)
    
    #plot fit, if wanted
       #plot(ln.avg.mass ~ D.elapse) #CHECK W/CAM REGARDING "[-9]"
       #abline(growth.lm)
    #   
    #Global Instaneaous growth rate (igr)
    #igr = as.numeric(coef(growth.lm)[2]);print(igr) 
    #intercept = coef(growth.lm)[1]
##### End old code section #####    
    
    #interval IGR
    #browser()
    IGR.tab = data.frame(nrow  = length(ln.avg.mass), ncol = 2 )
     for(i in 1:length(dates)){
       IGR.tab[i,1] = as.numeric(ln.avg.mass[i])
       IGR.tab[i,2] = as.numeric(ln.avg.mass[i+1])
       IGR.tab[i,3] = as.numeric(dates[i+1] - dates[i])
       IGR.tab[i,4] = abs(as.numeric(ln.avg.mass[i+1])-as.numeric(ln.avg.mass[i]))/as.numeric(IGR.tab[i,3])
     }
    IGR.tab = as.data.frame(IGR.tab[-nrow(IGR.tab),])
    colnames(IGR.tab) = c("ln.start.mass", "ln.end.mass", "d", "IGRint")
    rownames(IGR.tab) = NULL
   
    ## 4) Estimate P  
    P.tab = data.frame(ncol = 4, nrow = length(dates))
    Biomass = as.numeric(Biomass.tab[,"Biomass"])
    for (i in 1:length(dates)){
      P.tab[i,1] = as.character(dates[i])
      P.tab[i,2] = as.character(dates[i +1])
      P.tab[i,3] = as.numeric(dates[i+1] - dates[i])
      P.tab[i,4] = as.numeric(mean(Biomass[i:(i+1)])*as.numeric(IGR.tab[i,4])*as.numeric(P.tab[i,3]))
      }
    
    #clean up P.tab
    #P.tab = data.frame(P.tab)
    P.tab = as.data.frame(P.tab[-nrow(P.tab),])
    colnames(P.tab) = c("start_date", "end_date", "d", "Pint") 
    rownames(P.tab) = NULL
    
   #Calculate P and P/B (P_B), and Pint
    P = sum(P.tab[,"Pint"])
    P_B = P/mean(Biomass)
    
    #Calculate N and B
    N = mean(as.numeric(Biomass.tab[,"Abundance"]))
    B = mean(Biomass)
    Sums = c(stream, as.character(taxon$TAXON), cohort, N, B, P, P_B, igr)
    #print(Sums)
    MCsum[n,] = Sums 
    
    #need to automate this to length of dates
    Biomass.tab = data.frame(Biomass.tab)
    colnames(Biomass.tab)[c(1,(dim(Biomass.tab)[2]-3):(dim(Biomass.tab)[2]-2))] = c("start_date", "MASS_START", "BTOT_START")
    B.merge = merge(P.tab, Biomass.tab[,c(1,(dim(Biomass.tab)[2]-3):(dim(Biomass.tab)[2]-2))], by = "start_date")
    colnames(Biomass.tab)[c(1,(dim(Biomass.tab)[2]-3):(dim(Biomass.tab)[2]-2))] = c("end_date", "MASS_END", "BTOT_END")
    B.merge = merge(B.merge, Biomass.tab[,c(1,(dim(Biomass.tab)[2]-3):(dim(Biomass.tab)[2]-2))], by = "end_date")

    colnames(Biomass.tab)[c(1:2)] = c("start_date", "N_START")
    N.merge = merge(B.merge, Biomass.tab[,c(1:2)], by = "start_date")
    colnames(Biomass.tab)[c(1:2)] = c("end_date", "N_END")
    N.merge = merge(N.merge, Biomass.tab[,c(1:2)], by = "end_date")
    colnames(Biomass.tab)[c(1)] = c("start_date")
    N.merge = merge(N.merge, Biomass.tab[,c(1,dim(Biomass.tab)[2]-1,dim(Biomass.tab)[2])], by = "start_date")
  #browser()
    #print(N.merge$Pint)
    #print(N.merge$BTOT_START)
    #print(N.merge$BTOT_END)
 
    N.merge$P_Bint = apply(N.merge, 1, function(x) as.numeric(x['Pint'])/((as.numeric(x['BTOT_START']) + as.numeric(x['BTOT_END']))/2))
    N.merge$igr = IGR.tab[,4]

     ## 5) Wrap up output
    #print(as.matrix(N.merge))
    #browser()
    #print(rep(stream, nrow(N.merge))); print(rep(as.character(taxon$TAXON)));print(rep(cohort),nrow(N.merge))
    #print(N.merge)
    Pars = cbind(rep(stream, length(nrow(N.merge))), rep(as.character(taxon$TAXON), length(nrow(N.merge))),
                 rep(cohort, length(nrow(N.merge))), as.matrix(N.merge))
    #Pars = data.frame(stream,as.character(taxon$TAXON),cohort,N.merge)
    #print(Pars)
    #print(unlist(Pars))
    #MCsims[n,] = unlist(Pars)
    
    MCsims = rbind(MCsims, Pars)
  }#need to add MCsims clean up outside of this
  #clean up MCsum and MCsims
  rownames(MCsum) = NULL
  colnames(MCsum) = c("SITE", "TAXON", "COHORT", "N", "B", "P", "P_B", "IGR")
  rownames(MCsims) = NULL
  colnames(MCsims) =  c("SITE","TAXON","COHORT","END_DATE","START_DATE","DAYS", 
                     "P.INT", "IND_MASS_START", "BTOT_START",
                    "IND_MASS_END", "BTOT_END", "N_START","N_END","START_SIZE", "END_SIZE",
                    "P_B", "IGRint")
  #colnames(MCsims) =  c("SITE","TAXON","COHORT","END_DATE","START_DATE","DAYS", 
     #                   "P.INT", "IND_MASS_START", "MIN_SIZE", "BTOT_START",
      #                  "IND_MASS_END", "MAX_SIZE",
       #                 "BTOT_END", "N_START","N_END","P_B", "IGRint")
  
  #print(summary(MCsims))
  #print(MCsims)
  #plot summary
  #par(mfrow = c(length(Sums),1), mar = rep(0.7,4))
  #par(mfrow = c(5,1), mar = rep(1,4))
  #a= hist(MCsum[,"N"])
  #text(x= 0.95*(max(a$breaks)), y = 0.9*(max(a$counts)), paste("nSims = ", nBoot ))
  #hist(MCsum[,"B"])
  #hist(MCsum[,"P"])
  #hist(MCsum[,"igr"])
  #hist(MCsum[,"P_B"])
  #browser()
  #print(file.name.out)
  #print(paste("SUMMARY_",file.name.out, sep = ""))
  write.csv(MCsims, file = file.name.out, row.names = F)
  write.csv(MCsum, file = paste("SUMMARY_",gsub(x = a, '.csv', replacement = ""),"_Prod.csv", sep = ""), row.names = F)
  }
}


