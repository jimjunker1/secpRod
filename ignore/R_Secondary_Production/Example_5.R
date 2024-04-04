#This code calculates production for RM30, year2.5 (9/15/2007-9/14/2008)


#Set the working directory:
  #All calls are relative to this directory on your machine
  #Either use the 'setwd' command below (and insert the appropriate directory)
  #or use 'Change dir...' from the 'File' drop-down menu and browse to the appropriate directory
  setwd("C:/Users/Jim/Documents/Projects/General Files/R Secondary Production-GC Example")


#Libraries:
  #These packages must be installed on your machine:
  #  'chron'
  #Load the package that can handle dates:
  library(chron)


#Define functions stored in separate scripts:
  source("Scripts/get.dates_function.txt")
  source("Scripts/igr.prod_function.txt")
  source("Scripts/sf.prod_function.txt")
  source("Scripts/pb.prod_function.txt")
  source("Scripts/wrapper.site.yr_function.txt")
  source("Scripts/hab.weight_function.txt")
  source("Scripts/sampinfo.site.yr_function.txt")


#Import raw data:
  #Data is a matrix of abundance (number m^-2) for each size class (length in mm, columns) in each sample (rows)
  #Five additional columns preceed the data to denote SITE, SAMPLE, DATE, HABITAT, & TAXON
  #(This example dataset is for river mile 30 (all taxa, habitats, & years), and also has all Chironomidae data from river mile 0 in year 3 appended to the end)
  EX <- read.table("Data/RM_30_EXAMPLE.txt", header=T, sep="\t", quote="", strip.white=TRUE)

  #Extract years, months, and days, and re-code them as numeric variables:
  year <- as.numeric(as.character(years(chron(dates = as.character(EX$DATE)))))
  month <- as.numeric(months(chron(dates = as.character(EX$DATE))))
  day <- as.numeric(days(chron(dates = as.character(EX$DATE))))

  #Combine year, month, and day into a single julian date variable starting with January 1, 2006 (i.e.- Jan 1, 2006 = day1):
  JULIAN <- julian(month, day, year, origin=c(month = 12, day = 31, year = 2005))

  #Insert the julian date variable into the dataframe:
  EX <- data.frame(EX[,1:2], JULIAN, EX[,3:(dim(EX)[2])])


#Import temperature data:
  #Average temperatures (degrees C) for each of the 4 sampling intervals for RM30, year2.5
  #DON'T HAVE REAL TEMP DATA FOR RM30, SO I MADE IT UP BY SUBSETTING RM0 DATA:
  temp <- scan("Data/Chiro_temp.txt")
  temp <- temp[1:4]


#Import taxonomic info:
  #L-M parameters, method of production, growth paramters, cpi's, p/b's, etc. for each taxon
  tax <- read.table("Taxon Info/taxa_info_for_R.txt", header=T, sep="\t", quote="", strip.white=TRUE)



#See all the sampling dates for the specified site and habitat:
  get.dates(DATA=EX, site="RM30", habitat="COBBLE")



#See all the sampling dates bounded by the specified range for the specified site and habitat:
  get.dates(DATA=EX, site="RM30", habitat="COBBLE", first.date="9/15/2007", last.date="9/14/2008")



#Run the igr.prod function on the selected data:
  igr.out <- igr.prod(	DATA=EX,
			site="RM30", first.date="9/15/2007", last.date="9/14/2008",
			habitat="COBBLE", taxon="Chironomid",
			TEMP=temp, wrap=F,
			LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
			g.a=0.051, g.b=0.068, g.c=0.006, g.d=NA, min.growth=0, temp.corr=c(1,1,1,1),
			boot.num=100)

  names(igr.out)			#names of the objects in the output list
  igr.out[c(1,2,3,7,8,9,11,12)]		#looking only at the short ones in this list

  hist(igr.out$Pboots, xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="")
    abline(v=mean(igr.out$Pboots), lwd=3, col="red")
    abline(v=median(igr.out$Pboots), lwd=3, col="blue")
    abline(v=quantile(igr.out$Pboots, c(0.025, 0.975)), lwd=3, col="green")
    
    rm(igr.out)
    
    
#Look at the scaling of estimates with increasing boot size.
    
#Run the igr.prod function on the selected data:
    

    igr.out2 <- igr.prod(	DATA=EX,
                         site="RM30", first.date="9/15/2007", last.date="9/14/2008",
                         habitat="COBBLE", taxon="Chironomid",
                         TEMP=temp, wrap=F,
                         LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
                         g.a=0.051, g.b=0.068, g.c=0.006, g.d=NA, min.growth=0, temp.corr=c(1,1,1,1),
                         boot.num=200)
    
    names(igr.out2)			#names of the objects in the output list
    igr.out2[c(1,2,3,7,8,9,11,12)]		#looking only at the short ones in this list
    
    hist(igr.out2$Pboots, xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="")
    abline(v=mean(igr.out2$Pboots), lwd=3, col="red")
    abline(v=median(igr.out2$Pboots), lwd=3, col="blue")
    abline(v=quantile(igr.out2$Pboots, c(0.025, 0.975)), lwd=3, col="green")
    
    rm(igr.out2)
##500    

    igr.out5 <- igr.prod(	DATA=EX,
                          site="RM30", first.date="9/15/2007", last.date="9/14/2008",
                          habitat="COBBLE", taxon="Chironomid",
                          TEMP=temp, wrap=F,
                          LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
                          g.a=0.051, g.b=0.068, g.c=0.006, g.d=NA, min.growth=0, temp.corr=c(1,1,1,1),
                          boot.num=500)
    
    names(igr.out5)			#names of the objects in the output list
    igr.out5[c(1,2,3,7,8,9,11,12)]#looking only at the short ones in this list
    gc(igr.out5)
    rm(igr.out5)

#1000

    igr.out10 <- igr.prod(	DATA=EX,
                          site="RM30", first.date="9/15/2007", last.date="9/14/2008",
                          habitat="COBBLE", taxon="Chironomid",
                          TEMP=temp, wrap=F,
                          LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
                          g.a=0.051, g.b=0.068, g.c=0.006, g.d=NA, min.growth=0, temp.corr=c(1,1,1,1),
                          boot.num=1000)
    
    names(igr.out10)			#names of the objects in the output list
    igr.out10[c(1,2,3,7,8,9,11,12)]		#looking only at the short ones in this list
    gc(igr.out10)
    rm(igr.out10)
#10000
        igr.out100 <- igr.prod(	DATA=EX,
                           site="RM30", first.date="9/15/2007", last.date="9/14/2008",
                           habitat="COBBLE", taxon="Chironomid",
                           TEMP=temp, wrap=F,
                           LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
                           g.a=0.051, g.b=0.068, g.c=0.006, g.d=NA, min.growth=0, temp.corr=c(1,1,1,1),
                           boot.num=10000)
    
    names(igr.out100)			#names of the objects in the output list
    igr.out100[c(1,2,3,7,8,9,11,12)]		#looking only at the short ones in this list
    rm(igr.out100)
  gc(igr.out100)
#Run the sf.prod function on the selected data:
  sf.out <- sf.prod(	DATA=EX,
			site="RM30", first.date="9/15/2007", last.date="9/14/2008",
			habitat="COBBLE", taxon="Simuliid",
			LM.a=0.004, LM.b=2.807, LM.p.ash=0,
			num.size.classes=11, min.cpi=300, max.cpi=365, temp.corr=1,
			boot.num=100)

  names(sf.out)				#names of the objects in the output list
  sf.out[c(1,3,4,8,9,10,11)]		#looking only at the short ones in this list

  hist(sf.out$Pboots, xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="")
    abline(v=mean(sf.out$Pboots), lwd=3, col="red")
    abline(v=median(sf.out$Pboots), lwd=3, col="blue")
    abline(v=quantile(sf.out$Pboots, c(0.025, 0.975)), lwd=3, col="green")

    
## Run chironomids with sf method
    
sf.outC <- sf.prod(	DATA=EX,
                       site="RM30", first.date="9/15/2007", last.date="9/14/2008",
                       habitat="COBBLE", taxon="Chironomid",
                       LM.a=0.0006, LM.b=2.77, LM.p.ash=0,
                       num.size.classes=14, min.cpi=335, max.cpi=365, temp.corr=1,
                       boot.num=500)
    
    names(sf.outC)				#names of the objects in the output list
    sf.outC[c(1,3,4,8,9,10,11)]		#looking only at the short ones in this list
    
    hist(sf.outC$Pboots, xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="")
    abline(v=mean(sf.outC$Pboots), lwd=3, col="red")
    abline(v=median(sf.outC$Pboots), lwd=3, col="blue")
    abline(v=quantile(sf.outC$Pboots, c(0.025, 0.975)), lwd=3, col="green")


#Run the pb.prod function on the selected data:
  pb.out <- pb.prod(	DATA=EX,
			site="RM30", first.date="9/15/2007", last.date="9/14/2008",
			habitat="COBBLE", taxon="Mite",
			LM.a=0.00266, LM.b=1, LM.p.ash=0,
			p.b=5, temp.corr=1,
			boot.num=100)

  names(pb.out)				#names of the objects in the output list
  pb.out[c(1,2,3,4,8,9,10,11)]		#looking only at the short ones in this list

  hist(pb.out$Pboots, xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="")
    abline(v=mean(pb.out$Pboots), lwd=3, col="red")
    abline(v=median(pb.out$Pboots), lwd=3, col="blue")
    abline(v=quantile(pb.out$Pboots, c(0.025, 0.975)), lwd=3, col="green")



#Run the wrapper function on the selected data:
  rm30y2.5.out <- wrapper.site.yr(	DATA=EX,
					site="RM30", first.date="9/15/2007", last.date="9/14/2008",
					habitat="ALL", TAXA=tax,
					TEMP.COB=temp, TEMP.DEP=temp, TEMP.TAL=temp, wrap=F,
					temp.corr.igr.cob=c(1,1,1,1), temp.corr.igr.dep=c(1,1,1,1), temp.corr.igr.tal=c(1,1,1,1),
					temp.corr.sfpb=1, boot.num=50)

  names(rm30y2.5.out)				#names of the objects in the output list

  #Taking a closer look at only the production results:

  rm30y2.5.out$Pboots.cob[1:10,]		#looking at only the first 10 bootstrap estimates
  rm30y2.5.out$Pboots.dep[1:10,]		#looking at only the first 10 bootstrap estimates
  rm30y2.5.out$Pboots.tal[1:10,]		#looking at only the first 10 bootstrap estimates

  rm30y2.5.out$Pboots.cob[[2]]			#reference the output results like this...
  rm30y2.5.out$Pboots.cob$Chironomid		#or like this

  mean(rm30y2.5.out$Pboots.cob[[2]])					#Summary bootstrap estimates of Chironomids for COBBLE
  quantile(rm30y2.5.out$Pboots.cob[[2]], c(0.025, 0.5, 0.975))
  mean(rm30y2.5.out$Pboots.dep[[2]])					#Summary bootstrap estimates of Chironomids for DEPOSITIONAL
  quantile(rm30y2.5.out$Pboots.dep[[2]], c(0.025, 0.5, 0.975))
  mean(rm30y2.5.out$Pboots.tal[[2]])					#Summary bootstrap estimates of Chironomids for TALUS.CLIFF
  quantile(rm30y2.5.out$Pboots.tal[[2]], c(0.025, 0.5, 0.975))

  #Graph bootstrap estimates for this taxon in all 3 habitats:
  par(mfrow=c(3,1))
  hist(rm30y2.5.out$Pboots.cob[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="COBBLE")
    abline(v=mean(rm30y2.5.out$Pboots.cob[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.cob[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.cob[[2]], c(0.025, 0.975)), lwd=3, col="green")
  hist(rm30y2.5.out$Pboots.dep[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="DEPOSITIONAL")
    abline(v=mean(rm30y2.5.out$Pboots.dep[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.dep[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.dep[[2]], c(0.025, 0.975)), lwd=3, col="green")
  hist(rm30y2.5.out$Pboots.tal[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="TALUS.CLIFF")
    abline(v=mean(rm30y2.5.out$Pboots.tal[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.tal[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.tal[[2]], c(0.025, 0.975)), lwd=3, col="green")
  par(mfrow=c(1,1))

  #Graph bootstrap estimates for this taxon in all 3 habitats (axes same scale):
  par(mfrow=c(3,1))
  hist(rm30y2.5.out$Pboots.cob[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="COBBLE", xlim=c(0,max(c(rm30y2.5.out$Pboots.cob[[2]], rm30y2.5.out$Pboots.dep[[2]], rm30y2.5.out$Pboots.tal[[2]]))))
    abline(v=mean(rm30y2.5.out$Pboots.cob[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.cob[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.cob[[2]], c(0.025, 0.975)), lwd=3, col="green")
  hist(rm30y2.5.out$Pboots.dep[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="DEPOSITIONAL", xlim=c(0,max(c(rm30y2.5.out$Pboots.cob[[2]], rm30y2.5.out$Pboots.dep[[2]], rm30y2.5.out$Pboots.tal[[2]]))))
    abline(v=mean(rm30y2.5.out$Pboots.dep[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.dep[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.dep[[2]], c(0.025, 0.975)), lwd=3, col="green")
  hist(rm30y2.5.out$Pboots.tal[[2]], xlab=expression(paste("Annual Production (mg AFDM m"^-2, " yr"^-1, ")", sep="")), main="TALUS.CLIFF", xlim=c(0,max(c(rm30y2.5.out$Pboots.cob[[2]], rm30y2.5.out$Pboots.dep[[2]], rm30y2.5.out$Pboots.tal[[2]]))))
    abline(v=mean(rm30y2.5.out$Pboots.tal[[2]]), lwd=3, col="red")
    abline(v=median(rm30y2.5.out$Pboots.tal[[2]]), lwd=3, col="blue")
    abline(v=quantile(rm30y2.5.out$Pboots.tal[[2]], c(0.025, 0.975)), lwd=3, col="green")
  par(mfrow=c(1,1))



#Calculate habitat-weighted production for all taxa (habitat weights are made-up):
hab.weight(dat.cob=rm30y2.5.out$Pboots.cob, dat.dep=rm30y2.5.out$Pboots.dep, dat.tal=rm30y2.5.out$Pboots.tal, wt.cob=0.4, wt.dep=0.4, wt.tal=0.2)



#Calculate habitat-weighted production for all taxa (habitat weights are made-up):
sampinfo.site.yr(DATA=EX, site="RM30", first.date="9/15/2007", last.date="9/14/2008", habitat="ALL", TAXA=tax)



