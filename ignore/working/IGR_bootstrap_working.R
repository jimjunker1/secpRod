#2oP Bootstrap, Author: Cam MacKenzie; Modified by: Jim Junker
#cmackenzie@atwaterresources.com
#james.junker1@gmail.com
#last modified: Jan 2018
#Coded from 2oP macro 
#========================================================
# Input data 
#==============================================================================
# Data (.csv) MUST be in the following format, in particular the date (i.e. %m/%d/%y), 
# but number of columns before size class data should also be the same
# (#/S = #/sample or density (#/m-2):
# DATE      SAMPLE  High taxon  Low taxon #/S   1 2  3  4	5	6.....
# 6/18/2012	-0.42	  Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/18/2012	-0.42	  Baetidae	  Acen	    44	  44.44444	0	0	0	0	0
# 6/18/2012	-0.17	  Baetidae	  Acen	    0 	  0	0	0	0	0	0
# 6/18/2012	0	      Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/18/2012	-0.17	  Baetidae	  Acen	    44	  44.44444	0	0	0	0	0
# 6/18/2012	0	      Baetidae	  Acen	    0	    0	0	0	0	0	0
# 6/27/2012	-0.17	  Baetidae	  Acen	    1022	577.77772	444.4444	0	0	0	0
# 6/27/2012	-0.38	  Baetidae	  Acen	    311	  266.66664	44.44444	0	0	0	0
# 6/27/2012	-0.42	  Baetidae	  Acen	    889	  799.99992	88.88888	0	0	0	0
#...

# Output format
#===============================================================================
# Data (.csv) will be written with the following columns
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
#Be sure individual Mass values are appropriate for species
#  size classes used
#Format Date
#Set appropriate # of Production Intervals

#need to update so output the DATE column for igr estimate

graphics.off()
#rm(list=ls())
# set.seed(100) # sets randomization seed for bootstrapping, turn on if you want the same results every time you run the code 
#setwd("~/Rdata") #Set working directory

IGR_file_boot <- function(path, file.names,tax.info, nBoot){
  library(data.table)
  for(a in file.names){
  allData <- read.table(file = paste(path,a, sep = ""),header=TRUE, sep = "\t",check.names = F,stringsAsFactors=FALSE)
  file.name.out = paste(path,"/IGR_boot_output/",gsub(x = a, '.txt', replacement = ""),"_prod.csv", sep = "") 
  # Define samples to crop (due to outliers) for igr calculation  
  # enter 0 (don't crop) or crop samples as follows c(1,9) or 9;   
  # removes first and ninth sample or ninth sample only, respectively
  # for igr calculation 
  crop.dates = 0
  # Define number of Monte Carlo simulations to run
  nBoot= nBoot
  
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
source("./Analysis/functions/bootstrap_function.R")  
  bootsData = boots(allData)
  
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
    
    # loop to determine maximum size class, last column with data
    for (j in length(sum.classes):1){
      if (sum.classes[j]>0){
        no.classes = j
        break
      }
    }
    #can do this easier to speed up code in the future
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
      
      for (j in length(sum.classes):1){
        if (sum.classes[j]>0){
          classes = j
          break
        }
      }
      END_SIZE = as.numeric(colnames(size.classes[as.numeric(j)]))
      NonNAindex = which(!is.na(sum.classes))
      START_SIZE = as.numeric(colnames(size.classes[min(NonNAindex)]))
      for (k in 1:no.classes){
        #calculate mean of bootstrapped size class for date i, rounded as in macro
        N.tab[i, k+1] = as.numeric(sum(date.sub[i,k+5]))/length(date.sub[i,1]) #values are rounded as in macros
      }
      N.tab[i, no.classes+2] = START_SIZE
      N.tab[i, no.classes+3] = END_SIZE
    }
    rownames(N.tab) = NULL
    colnames(N.tab) = c("Date", 1:no.classes, "START_SIZE", "END_SIZE")
    N.tab[,1] = as.Date(N.tab[,1]) 
    ## 2) Calculate Biomass (B) table
   #need to quite dealin g with this vector of varying types crap
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
    ln.avg.mass = log(as.numeric(as.character(Biomass.tab[,"avg mass"])))
    
    #Crop dates for appropriate IGR model (specified by crop.dates above)
    if (crop.dates[1] != 0){
      D.elapse = D.elapse[-crop.dates]
      ln.avg.mass = ln.avg.mass[-crop.dates]
    }
    
    # Run linear model
    growth.lm = lm(ln.avg.mass ~ D.elapse)
    
    #summary.lm(growth.lm)
  #browser()
    # Diagnostic plots, if wanted
    #plot(growth.lm)
    
    #plot fit, if wanted
       #plot(ln.avg.mass ~ D.elapse) #CHECK W/CAM REGARDING "[-9]"
       #abline(growth.lm)
    #   
    #Global Instaneaous growth rate (igr)
    igr = as.numeric(coef(growth.lm)[2]) 
    intercept = coef(growth.lm)[1]
    
    #interval IGR
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
    MCsum[n,] = Sums 
    browser()
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
    colnames(Biomass.tab)[c(1:2)] = c("start_date", "N_START")
    N.merge = merge(N.merge, Biomass.tab[,c(1:2,dim(Biomass.tab)[2]-1,dim(Biomass.tab)[2])], by = "start_date")
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


######
#place to play with code before inserting

allData = read.csv(file = "./Output/stream_cohorts/Hver_Limnophora_A.csv", check.names = F)
i = "Hver_Limnophora_A.csv"
class(i)
grep(pattern = ".csv", x = i, value = T)


substr(names(allData)[6:(dim(allData)[2])], 1,6)
substr("abcdefg", 2,4)
########
#clean up MCsims 
rownames(MCsims) = NULL
colnames(MCsims) =  c("N","B","P","igr","P_B", "Pint.1", "Pint.2", "Pint.3", "Bint.1","Bint.2", "Bint.3", "Bint.4", "Nint.1","Nint.2", "Nint.3", "Nint.4", "size.1", "size.2", "size.3", "size.4")

print(MCsims)

#print summary
print(summary(MCsims))

#plot summary
par(mfrow = c(length(Pars),1), mar = rep(0.7,4))
par(mfrow = c(5,1), mar = rep(1,4))
a= hist(MCsims[,"N"])
text(x= 0.95*(max(a$breaks)), y = 0.9*(max(a$counts)), paste("nSims = ", nBoot ))
hist(MCsims[,"B"])
hist(MCsims[,"P"])
hist(MCsims[,"igr"])
hist(MCsims[,"P_B"])


##write to a file
#write.csv(MCsims, file = file.name.out )


int1 <- D.elapse[2] - D.elapse[1]#days in interval
int2 <- D.elapse[3] - D.elapse[2]
int3 <- D.elapse[4] - D.elapse[3]


Pd1 <- MCsims[,"Pint.1"]/int1#interval production/days in interval (P/d)
Pd2 <- MCsims[,"Pint.2"]/int2
Pd3 <- MCsims[,"Pint.3"]/int3

Pd = data.frame(Pd1 = Pd1, Pd2 = Pd2, Pd3 = Pd3)


Bavg1 <- (MCsims[, "Bint.1"] + MCsims[, "Bint.2"]) / 2#Average biomass between Sample date 1 & 2
Bavg2 <- (MCsims[, "Bint.2"] + MCsims[, "Bint.3"]) / 2
Bavg3 <- (MCsims[, "Bint.3"] + MCsims[, "Bint.4"]) / 2


Pd.B1 <- Pd1/Bavg1#interval production/mean interval biomass or daily growth rate g/g*d
Pd.B2 <- Pd2/Bavg2
Pd.B3 <- Pd3/Bavg3

Pd.B = data.frame(Pd.B1 = Pd.B1, Pd.B2 = Pd.B2, Pd.B3 = Pd.B3)#dataframe of booted interval P/B for all intervals

#Function for 95% CIs
CI <- function(x)
{
quantile(x, c(0.025, 0.975))
}
#Calculate the column CIs
CI95.1 <- apply(MCsims, 2, CI)
CI95.2 <- apply(Pd, 2, CI)
CI95.3 <- apply(Pd.B, 2, CI)

#Calculate the column means
boot.mean.1 <- apply(MCsims, 2, mean)
boot.mean.2 <- apply(Pd, 2, mean)
boot.mean.3 <- apply(Pd.B, 2, mean)

total1 <- rbind(boot.mean.1, CI95.1)
total2 <- rbind(boot.mean.2, CI95.2)
total3 <- rbind(boot.mean.3, CI95.3)

q <- cbind(total1, total2)
r <- cbind(q, total3)

a <- t(r)
write.csv(a, file = file.name.out )

`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````


#mean(MCsims[,1])
#quantile(MCsims[,1], c(0.025, 0.975))

P.tab
d1 <- rep(0, each=1000)
d2 <- rep(33, each=1000)
df <- data.frame(days=c(d1, d2))
a <- MCsims[,17]
b<-MCsims[,18]
c <- data.frame(size=c(a, b))
plot(c$size~df$days)
reg <- lm(c$size~df$days)
abline(reg)
confint(reg)

pdf("sample.pdf", 7, 5)
plot(ln.avg.mass ~ D.elapse)
abline(growth.lm)
dev.off()



int1 <- D.elapse[2] - D.elapse[1]
int2 <- D.elapse[3] - D.elapse[2]
int3 <- D.elapse[4] - D.elapse[3]
int4 <- D.elapse[5] - D.elapse[4]
int5 <- D.elapse[6] - D.elapse[5]
int6 <- D.elapse[7] - D.elapse[6]
int7 <- D.elapse[8] - D.elapse[7]
int8 <- D.elapse[9] - D.elapse[8]
int9 <- D.elapse[10] - D.elapse[9]
int10 <- D.elapse[11] - D.elapse[10]

Pd1 <- MCsims[,"Pint.1"]/int1
Pd2 <- MCsims[,"Pint.2"]/int2
Pd3 <- MCsims[,"Pint.3"]/int3
Pd4 <- MCsims[,"Pint.4"]/int4
Pd5 <- MCsims[,"Pint.5"]/int5
Pd6 <- MCsims[,"Pint.6"]/int6
Pd7 <- MCsims[,"Pint.7"]/int7
Pd8 <- MCsims[,"Pint.8"]/int8
Pd9 <- MCsims[,"Pint.9"]/int9
Pd10 <- MCsims[,"Pint.10"]/int10
Pd = data.frame(Pd1 = Pd1, Pd2 = Pd2, Pd3 = Pd3, Pd4 = Pd4, Pd5 = Pd5, Pd6 = Pd6, Pd7 = Pd7, Pd8 = Pd8, Pd9 = Pd9, Pd10 = Pd10)


Bavg1 <- (MCsims[, "Bint.1"] + MCsims[, "Bint.2"]) / 2
Bavg2 <- (MCsims[, "Bint.2"] + MCsims[, "Bint.3"]) / 2
Bavg3 <- (MCsims[, "Bint.3"] + MCsims[, "Bint.4"]) / 2
Bavg4 <- (MCsims[, "Bint.4"] + MCsims[, "Bint.5"]) / 2
Bavg5 <- (MCsims[, "Bint.5"] + MCsims[, "Bint.6"]) / 2
Bavg6 <- (MCsims[, "Bint.6"] + MCsims[, "Bint.7"]) / 2
Bavg7 <- (MCsims[, "Bint.7"] + MCsims[, "Bint.8"]) / 2
Bavg8 <- (MCsims[, "Bint.8"] + MCsims[, "Bint.9"]) / 2
Bavg9 <- (MCsims[, "Bint.9"] + MCsims[, "Bint.10"]) / 2
Bavg10 <- (MCsims[, "Bint.10"] + MCsims[, "Bint.11"]) / 2

Pd.B1 <- Pd1/Bavg1
Pd.B2 <- Pd2/Bavg2
Pd.B3 <- Pd3/Bavg3
Pd.B4 <- Pd4/Bavg4
Pd.B5 <- Pd5/Bavg5
Pd.B6 <- Pd6/Bavg6
Pd.B7 <- Pd7/Bavg7
Pd.B8 <- Pd8/Bavg8
Pd.B9 <- Pd9/Bavg9
Pd.B10 <- Pd10/Bavg10
Pd.B = data.frame(Pd.B1 = Pd.B1, Pd.B2 = Pd.B2, Pd.B3 = Pd.B3, Pd.B4 = Pd.B4, Pd.B5 = Pd.B5, Pd.B6 = Pd.B6, Pd.B7 = Pd.B7, Pd.B8 = Pd.B8, Pd.B9 = Pd.B9, Pd.B10 = Pd.B10)
dailyP <- merge(Pd, Pd.B)



# for loop to calculate the number of days in each interval
days = NULL
for(n in 1:length(D.elapse) - 1) {
days[n] = D.elapse[n + 1] - D.elapse[n]
}


#putting the R functions for produciton here
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/get.dates_function.txt")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/hvigr.prod_function.R")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/sf.prod_function.txt")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/pb.prod_function.txt")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/wrapper.site.yr_function.R")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/hab.weight_function.txt")
source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/sampinfo.site.yr_function.txt")


tax = read.table(file = "./SFS code/hver_taxa_info_LISA_mod.txt", header = T, sep = "\t", quote = "", strip.white = T)
colnames(tax) = c("TAXON", 	"METHOD", 	"LM.a", "LM.b",	"LM.p.ash",	"g.a",	"g.b",	"g.c",	"g.d",	"min.cpi",	"max.cpi",	"num.size.classes",	"p.b",	"Growth.equation", 	"min.growth",	"notes")
#source("C:/Users/Jim/Documents/Projects/Talk/SFS 2017/SFS2017/R Secondary Production-GC Example/Scripts/wrapper.site.yr_function.txt")
set.seed(123)
hver.out = wrapper.site.yr(DATA = hver_bugs, site = "Hver", habitat = "COBBLE", TEMP.COB = hv_temp1, TEMP.DEP = hv_temp1, TEMP.TAL = hv_temp1, first.date = "08/02/11", last.date = "07/26/12",
                           TAXA = tax, temp.corr.igr.cob = c(1,1,1,1,1,1,1,1,1,1), temp.corr.igr.dep = c(1,1,1,1,1,1,1,1,1,1), temp.corr.igr.tal = c(1,1,1,1,1,1,1,1,1,1), wrap = T, boot.num = 50)

