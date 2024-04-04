
#Define IGR production bootstrap function:
  int.CED <- function(DATA, site, first.date, last.date, habitat, taxon, TEMP, wrap=F, LM.a, LM.b, LM.p.ash, g.a, g.b, g.c, g.d, min.growth=0.0001, temp.corr, boot.num){
    #Select all data (across all taxa) for the specified site, dates, & habitat (this allows us to get all possible dates & sample numbers for taxa not occuring in all samples on all dates):
    data1 <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == habitat,]
    data1 <- data1[!is.na(apply(data1[,7:(dim(data1)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
    row.names(data1) <- NULL							#Must reset row names to be sequential

    #Select only the subset of the data (site, dates, habitat, taxon) to calculate production for:
    data2 <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == habitat & DATA$TAXON == taxon,]
    data2 <- data2[!is.na(apply(data2[,7:(dim(data2)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
    row.names(data2) <- NULL							#Must reset row names to be sequential

    #Create a vector of sampling dates (in Julian units):
    t <- sort(union(data1$JULIAN, NULL))#check length of t here
    
    #Create a vector of durations (day) between successive sampling intervals:
    t.int <- diff(t)
    if (wrap==T){						#If the user wants to calculate a "wrap-around" interval for the last interval...
      t.int <- c(t.int, (365-sum(t.int)))			#Create an extra interval set to be of a duration so that all intervals sum to 365 days
    }
#Make sure that the user supplied the correct number of interval temperatures and temperature-correction factors:
    if (length(TEMP) == length(t.int) & length(temp.corr) == length(t.int)){

      #Create the appropriate samples of all "zeroes" for those dates and samples in which the taxon was not found:
      for (d in 1:length(t)){
        all.samps <- as.numeric(levels(factor(data1$SAMPLE[data1$JULIAN == t[d]])))
        act.samps <- as.numeric(levels(factor(data2$SAMPLE[data2$JULIAN == t[d]])))
        zeroes <- setdiff(all.samps, act.samps)
          if (length(all.samps) > length(act.samps)){
            for (z in 1:length(zeroes)){
              to.add <- data.frame(matrix(c(site, zeroes[z], t[d], NA, habitat, taxon, rep(0, (dim(data2)[2]-6))), 1, dim(data2)[2]))
              names(to.add) <- names(data2)
              if (nrow(data2) == 0){
                data2 <- to.add
              }
              else data2 <- rbind(data2, to.add)
            }
          }
      }
      for (i in c(2,3,7:dim(data2)[2])){
        data2[,i] <- as.numeric(as.character(data2[,i]))				#Must force the appropriate columns to be numeric
      }
print(taxon)#;browser()
#Create a dataframe of sizes:
mm <- as.numeric(names(data2)[7:(dim(data2)[2])])
if(taxon == "Sperchon glandulosus" | taxon == "Sperchon sp."){
  AFDM <- exp(LM.a+LM.b*log(mm))*(1-(LM.p.ash/100))
} else if(taxon == "Ostracoda"){
  AFDM <- 10^((LM.a +LM.b*log10(mm))*(1-(LM.p.ash/100)))
} else if(taxon == "Copepoda"){
  AFDM <- ((1-(LM.p.ash/100))*LM.a*(mm^LM.b))/1000
}else{
  AFDM <- (1-(LM.p.ash/100))*LM.a*(mm^LM.b)
}
sizes <- data.frame(mm, AFDM)

      #Calculate SAMPLE annual production for each size class:
      N <- matrix(0, length(t), length(sizes$mm))
      Bdatesampinfo <- Ndatesampinfo <- data.frame(matrix(0, length(t), 5))
      names(Bdatesampinfo) <- names(Ndatesampinfo) <- c("DATE","JULIAN","N","Mean","Stdev")
      for (d in 1:length(t)){
        rows <- data2[as.numeric(row.names(data2))[data2$JULIAN == t[d]],(7:dim(data2)[2])]
        N[d,] <- apply(rows, 2, mean)											#For each date, compute the mean abundance from the selected samples for each size class
        Ndatesampinfo[d,3:5] <- c(length(apply(rows, 1, sum)), mean(apply(rows, 1, sum)),sd(apply(rows, 1, sum)))	#For each date, also compute the total abundance (over all size classes) for each sample, and calculate the number of samples, mean, and standard deviation of the samples on that date
        rows.B <- t(t(rows)*sizes$AFDM)											#Same as rows (above) but biomass instead of abundance for each size class and sample
        Bdatesampinfo[d,3:5] <- c(length(apply(rows.B, 1, sum)), mean(apply(rows.B, 1, sum)),sd(apply(rows.B, 1, sum)))	#For each date, also compute the total biomass (over all size classes) for each sample, and calculate the number of samples, mean, and standard deviation of the samples on that date
      }
      Bdatesampinfo[,1:2] <- Ndatesampinfo[,1:2] <- get.dates(DATA=DATA, site=site, habitat=habitat, first.date=first.date, last.date=last.date)
      int.N <- (N[1:(length(t)-1),] + N[2:length(t),])/2		#Calculate interval abundance (take the mean of successive dates)
      int.B <- t(t(int.N)*sizes$AFDM)					#Calculate interval biomass (interval abundance * mass for each size class)
      if (wrap==T){							#If the user wants to calculate a "wrap-around" interval for the last interval...
        int.B <- rbind(int.B, (int.B[1,] + int.B[(dim(int.B)[1]),])/2)	#Add a row to represent interval biomass for an additional, final interval (use the mean interval biomass of the first and last date to get the biomass for the last interval)
        int.N <- rbind(int.N, (int.N[1,] + int.N[(dim(int.N)[1]),])/2)	#Add a row to represent interval abundance for an additional, final interval (use the mean interval abundance of the first and last date to get the abundance for the last interval)
      }
      int.CED <- matrix(0, length(t.int), length(sizes$mm))
      for (d in 1:length(t)){
          for (s in 1:length(sizes$mm)){
          for (f in 1:length(t.int)){
            int.CED[f,s] <- int.B[f,s]*((t.int[f]*1440)*(exp(2.84+(-0.25*log(sizes$AFDM[s]))+(-0.65/TEMP[d]))))/24.1
        } 
          }
      }
      int.CED[int.CED < 0] <- 0	
      # browser()
      #Set negative interval production values to zero
      CEDintsampest <- apply(int.CED, 1, sum)	#Sum all size classes to calculate annual production for each interval
      if (wrap==T){				#If the user wants to calculate a "wrap-around" interval for the last interval...
        int.B <- int.B[1:(dim(int.B)[1]-1),]	#Remove the row representing interval biomass for the additional, final "wrap-around" interval before exiting the current bootstrap step
        int.N <- int.N[1:(dim(int.N)[1]-1),]	#Remove the row representing interval abundance for the additional, final "wrap-around" interval before exiting the current bootstrap step
      }
      CEDsampest <- sum(CEDintsampest)		#Calculate SAMPLE annual production for all intervals combined

      #Generate bootstrap production vectors for each interval:
      CEDintboots <- matrix(NA, boot.num, length(t.int))
      Bintboots <- Nintboots <- matrix(NA, boot.num, (length(t)-1))
      Bdateboots <- Ndateboots <- matrix(NA, boot.num, length(t))
      for (b in 1:boot.num){
        N <- matrix(0, length(t), length(sizes$mm))
        for (d in 1:length(t)){
          num.samps <- length(as.numeric(row.names(data2))[data2$JULIAN == t[d]])				#Set the number of samples that exist for this date
          row.nums <- sample(as.numeric(row.names(data2))[data2$JULIAN == t[d]], size=num.samps, replace=T)	#For each date, select (with replacement) <num.samps> entire samples (where each sample selected includes the data for all size classes, and <num.samps> is the total number of samples that exist for that date)
          rows <- data2[row.nums,(7:dim(data2)[2])]
          N[d,] <- apply(rows, 2, mean)										#and then compute the mean abundance from the selected samples for each size class
        }
        int.N <- (N[1:(length(t)-1),] + N[2:length(t),])/2		#Calculate interval abundance (take the mean of successive dates)
        int.B <- t(t(int.N)*sizes$AFDM)					#Calculate interval biomass (interval abundance * mass for each size class)
        B <- t(t(N)*sizes$AFDM)						#Also calculate biomass (sample date abundance * mass for each size class)
        if (wrap==T){									#If the user wants to calculate a "wrap-around" interval for the last interval...
          int.B <- rbind(int.B, (int.B[1,] + int.B[(dim(int.B)[1]),])/2)		#Add a row to represent interval biomass for an additional, final interval (use the mean interval biomass of the first and last date to get the biomass for the last interval)
          int.N <- rbind(int.N, (int.N[1,] + int.N[(dim(int.N)[1]),])/2)		#Add a row to represent interval abundance for an additional, final interval (use the mean interval abundance of the first and last date to get the abundance for the last interval)
        }
        int.CED <- matrix(0, length(t.int), length(sizes$mm))
            for (d in 1:length(t)){
            for (s in 1:length(sizes$mm)){
              for (f in 1:length(t.int)){
                int.CED[f,s] <- int.B[f,s]*((t.int[f]*1440)*(exp(2.84+(-0.25*log(sizes$AFDM[s]))+(-0.65/TEMP[d]))))/24.1
              }
            }
          }
        int.CED[int.CED < 0] <- 0				#Set negative interval production values to zero
        CEDintboots[b,] <- apply(int.CED, 1, sum, na.rm = T)		#Sum all size classes to calculate annual production for each interval
        if (wrap==T){					#If the user wants to calculate a "wrap-around" interval for the last interval...
          int.B <- int.B[1:(dim(int.B)[1]-1),]		#Remove the row to representing interval biomass for the additional, final "wrap-around" interval before exiting the current bootstrap step
          int.N <- int.N[1:(dim(int.N)[1]-1),]		#Remove the row to representing interval abundance for the additional, final "wrap-around" interval before exiting the current bootstrap step
        }
        Bintboots[b,] <- apply(int.B, 1, sum, na.rm = T)		#Sum all size classes to calculate biomass for each interval
        Nintboots[b,] <- apply(int.N, 1, sum, na.rm = T)		#Sum all size classes to calculate abundance for each interval
        Bdateboots[b,] <- apply(B, 1, sum, na.rm = T)		#Calculate bootstrap matrix of biomass estimates on each sample date
        Ndateboots[b,] <- apply(N, 1, sum, na.rm = T)		#Calculate bootstrap matrix of abundance estimates on each sample date
      }
      CEDboots <- apply(CEDintboots, 1, sum, na.rm = T)			#Calculate bootstrap vector of annual production for all intervals combined
      Bboots <- apply(Bintboots, 1, mean, na.rm = T)			#Calculate bootstrap vector of annual biomass by taking the mean of all intervals
      Nboots <- apply(Nintboots, 1, mean, na.rm = T)			#Calculate bootstrap vector of annual abundance by taking the mean of all intervals
      list(totdays=sum(t.int), CEDsampest=CEDsampest, CEDbootest=c(Mean=mean(CEDboots), quantile(CEDboots, c(0.025, 0.5, 0.975))), CEDboots=CEDboots, Bboots=Bboots, Nboots=Nboots, julians=t, intdays=t.int, CEDintsampest=CEDintsampest, CEDintboots=CEDintboots, Bintboots = Bintboots, Nintboots = Nintboots, Bdatesampinfo=Bdatesampinfo, Ndatesampinfo=Ndatesampinfo, Bdateboots=Bdateboots, Ndateboots=Ndateboots)
    }

    #If the correct number of interval temperatures and temperature-correction factors was not supplied, print an error message:
    else if (length(TEMP) != length(temp.corr)){
      print("Error: the number of temperatures does not match the number of temperature-correction factors")
    }
    else if (length(TEMP) != length(t.int)){
      print("Error: the number of temperatures does not match the number of intervals")
    }
  }

