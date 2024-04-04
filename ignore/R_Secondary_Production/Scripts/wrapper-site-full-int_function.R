
#Define a wrapper function to run production calcs on all taxa and all habitats for a given site and year:
  wrapper.site.int <- function(DATA, site, first.date, last.date, habitat, TAXA, TEMP.COB, TEMP.DEP, TEMP.TAL, wrap=F, temp.corr.igr.cob, temp.corr.igr.dep, temp.corr.igr.tal, temp.corr.sfpb=1, boot.num, GROWTH_df,...){
    
    # full_lists = vector(mode = "list", length = 3)
    # names(full_lists) <- c('Nintboots.full','Bintboots.full','Pintboots.full')
    Nintboots.full <- Bintboots.full<- Pintboots.full<- list()
    
    if (habitat == "COBBLE" | habitat == "DEPOSITIONAL" | habitat =="TALUS.CLIFF"){
      data3 <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == habitat,]
      data3 <- data3[!is.na(apply(data3[,7:(dim(data3)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
      row.names(data3) <- NULL							#Reset row names to be sequential
      num.intervals <- length(levels(factor(data3$JULIAN)))-1			#Define the number of intervals
      if (wrap==T){
        num.intervals <- num.intervals + 1
      }
      num.intervals.cob <- num.intervals.dep <- num.intervals.tal <- "empty"	#Cannot define the number of intervals for all habitats when only one habitat is called
      taxa.levels <- levels(factor(data3$TAXON))				#Pull out the names of all taxa represented in the selected data
    }
    else if (habitat == "ALL"){
      data3 <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1],]
      data3 <- data3[!is.na(apply(data3[,7:(dim(data3)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
      row.names(data3) <- NULL							#Reset row names to be sequential
      num.intervals <- "empty"							#Cannot define the number of intervals for all habitats simultaneously
      taxa.levels <- levels(factor(data3$TAXON))				#Pull out the names of all taxa represented in the selected data
    }
    else print("Error: invalid habitat")
    
    #Make sure that the number of temperatures, temperature correction factors, and sampling intervals match:
    if (habitat == "COBBLE" & length(TEMP.COB) != length(temp.corr.igr.cob)){
      print("Error: the number of temperatures does not match the number of temperature-correction factors")
    }
    else if (habitat == "DEPOSITIONAL" & length(TEMP.DEP) != length(temp.corr.igr.dep)){
      print("Error: the number of temperatures does not match the number of temperature-correction factors")
    }
    else if (habitat == "TALUS.CLIFF" & length(TEMP.TAL) != length(temp.corr.igr.tal)){
      print("Error: the number of temperatures does not match the number of temperature-correction factors")
    }
    else if (habitat == "COBBLE" & num.intervals != length(TEMP.COB)){
      print("Error: the number of temperatures does not match the number of intervals. WTF")
      print(get.dates(DATA=DATA, site=site, habitat=habitat, first.date=first.date, last.date=last.date))
    }
    else if (habitat == "DEPOSITIONAL" & num.intervals != length(TEMP.DEP)){
      print("Error: the number of temperatures does not match the number of intervals")
      print(get.dates(DATA=DATA, site=site, habitat=habitat, first.date=first.date, last.date=last.date))
    }
    else if (habitat == "TALUS.CLIFF" & num.intervals != length(TEMP.TAL)){
      print("Error: the number of temperatures does not match the number of intervals")
      print(get.dates(DATA=DATA, site=site, habitat=habitat, first.date=first.date, last.date=last.date))
    }
    else if (habitat == "ALL" & (length(TEMP.COB) != length(temp.corr.igr.cob) | length(TEMP.DEP) != length(temp.corr.igr.dep) | length(TEMP.TAL) != length(temp.corr.igr.tal))){
      print("Error: the number of temperatures does not match the number of temperature-correction factors for at least one of the habitats")
    }
    else if (habitat =="ALL" & length(TEMP.COB) == length(temp.corr.igr.cob) & length(TEMP.DEP) == length(temp.corr.igr.dep) & length(TEMP.TAL) == length(temp.corr.igr.tal)){
      #Define the sampling dates & number of intervals for "COBBLE":
      data.cob <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == "COBBLE",]
      data.cob <- data.cob[!is.na(apply(data.cob[,7:(dim(data.cob)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
      row.names(data.cob) <- NULL								#Reset row names to be sequential
      num.intervals.cob <- length(levels(factor(data.cob$JULIAN)))-1				#Define the number of intervals
      if (wrap==T){
        num.intervals.cob <- num.intervals.cob + 1
      }
      #Define the sampling dates & number of intervals for "DEPOSITIONAL":
      data.dep <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == "DEPOSITIONAL",]
      data.dep <- data.dep[!is.na(apply(data.dep[,7:(dim(data.dep)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
      row.names(data.dep) <- NULL								#Reset row names to be sequential
      num.intervals.dep <- length(levels(factor(data.dep$JULIAN)))-1				#Define the number of intervals
      if (wrap==T){
        num.intervals.dep <- num.intervals.dep + 1
      }
      #Define the sampling dates & number of intervals for "TALUS.CLIFF":
      data.tal <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == "TALUS.CLIFF",]
      data.tal <- data.tal[!is.na(apply(data.tal[,7:(dim(data.tal)[2])], 1, sum)),]		#Remove rows with NA's (missing data)
      row.names(data.tal) <- NULL								#Reset row names to be sequential
      num.intervals.tal <- length(levels(factor(data.tal$JULIAN)))-1				#Define the number of intervals
      if (wrap==T){
        num.intervals.tal <- num.intervals.tal + 1
      }
    }
    if (habitat == "ALL" & (num.intervals.cob != length(TEMP.COB) | num.intervals.dep != length(TEMP.DEP) | num.intervals.tal != length(TEMP.TAL))){
      print("Error: the number of temperatures does not match the number of intervals for at least one of the habitats")
      print("Sample dates for COBBLE:")
      print(get.dates(DATA=DATA, site=site, habitat="COBBLE", first.date=first.date, last.date=last.date))
      print("Sample dates for DEPOSITIONAL:")
      print(get.dates(DATA=DATA, site=site, habitat="DEPOSITIONAL", first.date=first.date, last.date=last.date))
      print("Sample dates for TALUS.CLIFF:")
      print(get.dates(DATA=DATA, site=site, habitat="TALUS.CLIFF", first.date=first.date, last.date=last.date))
    }

    #If temps, temp correction factors, and number of intervals do not match exactly for each habitat, then do not run the "guts" of the function:
    if ((habitat == "COBBLE" & (num.intervals != length(TEMP.COB) | length(TEMP.COB) != length(temp.corr.igr.cob))) | 
      (habitat == "DEPOSITIONAL" & (num.intervals != length(TEMP.DEP) | length(TEMP.DEP) != length(temp.corr.igr.dep))) | 
      (habitat == "TALUS.CLIFF" & (num.intervals != length(TEMP.TAL) | length(TEMP.TAL) != length(temp.corr.igr.tal))) | 
      (habitat == "ALL" & (num.intervals.cob != length(TEMP.COB) | num.intervals.dep != length(TEMP.DEP) | num.intervals.tal != length(TEMP.TAL) | 
        length(TEMP.COB) != length(temp.corr.igr.cob) | length(TEMP.DEP) != length(temp.corr.igr.dep) | length(TEMP.TAL) != length(temp.corr.igr.tal)))){
	print("Error: length of temperature or temp.corr does not equal num.intervals")
    }

    #Otherwise, do run the "guts of the function:
    else{
    Nboots.tal <- Nboots.dep <- Nboots.cob <- Bboots.tal <- Bboots.dep <- Bboots.cob <- Pboots.tal <- Pboots.dep <- Pboots.cob <- data.frame(matrix(NA, nrow=boot.num, ncol=length(taxa.levels)))
    Pintboots.cob <- data.frame(matrix(NA, nrow = length(taxa.levels), ncol = num.intervals))
    Bintboots.cob <- Nintboots.cob <- data.frame(matrix(NA, nrow = length(taxa.levels), ncol = num.intervals-1))
    for (a in 1:length(taxa.levels)){
      if (!is.element(taxa.levels[a], TAXA$TAXON)){
        tax.info <- TAXA[1,]
        tax.info[1] <- taxa.levels[a]
        tax.info[2:dim(tax.info)[2]] <- NA
      }
      else {
        tax.info <- TAXA[TAXA$TAXON == taxa.levels[a],]
      }
      if (is.na(tax.info$METHOD)){
        print(paste("Error: no method of production is specified for this taxon: ", taxa.levels[a]))
        out <- out.cob <- out.dep <- out.tal <- data.frame(Pboots=rep(NA, boot.num), Bboots=rep(NA, boot.num), Nboots=rep(NA, boot.num))
        Nintboots.full[[a]] = NA;Bintboots.full[[a]]= NA;Pintboots.full[[a]] = NA;
        
      } else if (tax.info$METHOD == "igr"){#need to integrate the growth_df dataset in here and above.And just spit the full temp file inthere.
        if (habitat == "COBBLE"){
          out <- igr.int.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat=habitat,
                   taxon=as.character(tax.info$TAXON), TEMP=TEMP.COB, wrap=wrap,
                   LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                   g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                   temp.corr=temp.corr.igr.cob, boot.num=boot.num, GROWTH_df = GROWTH_df)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_", habitat, "_", taxa.levels[a], ".txt", sep=""))
          print(out)
          sink()
        }
        if (habitat == "DEPOSITIONAL"){
          out <- igr.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat=habitat,
                   taxon=as.character(tax.info$TAXON), TEMP=TEMP.DEP, wrap=wrap,
                   LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                   g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                   temp.corr=temp.corr.igr.dep, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_", habitat, "_", taxa.levels[a], ".txt", sep=""))
          print(out)
          sink()
        }
        if (habitat =="TALUS.CLIFF"){
          out <- igr.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat=habitat,
                   taxon=as.character(tax.info$TAXON), TEMP=TEMP.TAL, wrap=wrap,
                   LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                   g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                   temp.corr=temp.corr.igr.tal, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_", habitat, "_", taxa.levels[a], ".txt", sep=""))
          print(out)
          sink()
        }
        if (habitat == "ALL"){
          out.cob <- igr.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="COBBLE",
                       taxon=as.character(tax.info$TAXON), TEMP=TEMP.COB, wrap=wrap,
                       LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                       temp.corr=temp.corr.igr.cob, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_COBBLE_", taxa.levels[a], ".txt", sep=""))
          print(out.cob)
          sink()
          out.dep <- igr.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="DEPOSITIONAL",
                       taxon=as.character(tax.info$TAXON), TEMP=TEMP.DEP, wrap=wrap,
                       LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                       temp.corr=temp.corr.igr.dep, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_DEPOSITIONAL_", taxa.levels[a], ".txt", sep=""))
          print(out.dep)
          sink()
          out.tal <- igr.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="TALUS.CLIFF",
                       taxon=as.character(tax.info$TAXON), TEMP=TEMP.TAL, wrap=wrap,
                       LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       g.a=tax.info$g.a, g.b=tax.info$g.b, g.c=tax.info$g.c, g.d=tax.info$g.d, min.growth=tax.info$min.growth,
                       temp.corr=temp.corr.igr.tal, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_TALUS.CLIFF_", taxa.levels[a], ".txt", sep=""))
          print(out.tal)
          sink()
        }
      }
      else if (tax.info$METHOD == "sf"){
        if (habitat == "COBBLE" | habitat == "DEPOSITIONAL" | habitat =="TALUS.CLIFF"){
          out <- sf.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat=habitat,
                   taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                   num.size.classes=tax.info$num.size.classes, min.cpi=tax.info$min.cpi, max.cpi=tax.info$max.cpi,
                   temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_", habitat, "_", taxa.levels[a], ".txt", sep=""))
          print(out)
          sink()
        }
        if (habitat == "ALL"){
          out.cob <- sf.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="COBBLE",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       num.size.classes=tax.info$num.size.classes, min.cpi=tax.info$min.cpi, max.cpi=tax.info$max.cpi,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_COBBLE_", taxa.levels[a], ".txt", sep=""))
          print(out.cob)
          sink()
          out.dep <- sf.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="DEPOSITIONAL",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       num.size.classes=tax.info$num.size.classes, min.cpi=tax.info$min.cpi, max.cpi=tax.info$max.cpi,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_DEPOSITIONAL_", taxa.levels[a], ".txt", sep=""))
          print(out.dep)
          sink()
          out.tal <- sf.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="TALUS.CLIFF",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash,
                       num.size.classes=tax.info$num.size.classes, min.cpi=tax.info$min.cpi, max.cpi=tax.info$max.cpi,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_TALUS.CLIFF_", taxa.levels[a], ".txt", sep=""))
          print(out.tal)
          sink()
        }
      }
      else if (tax.info$METHOD == "pb"){
        if (habitat == "COBBLE" | habitat == "DEPOSITIONAL" | habitat =="TALUS.CLIFF"){
          out <- pb.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat=habitat,
                   taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash, p.b=tax.info$p.b,
                   temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_", habitat, "_", taxa.levels[a], ".txt", sep=""))
          print(out)
          sink()
        }
        if (habitat == "ALL"){
          out.cob <- pb.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="COBBLE",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash, p.b=tax.info$p.b,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_COBBLE_", taxa.levels[a], ".txt", sep=""))
          print(out.cob)
          sink()
          out.dep <- pb.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="DEPOSITIONAL",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash, p.b=tax.info$p.b,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_DEPOSITIONAL_", taxa.levels[a], ".txt", sep=""))
          print(out.dep)
          sink()
          out.tal <- pb.prod(DATA=DATA, site=site, first.date=first.date, last.date=last.date, habitat="TALUS.CLIFF",
                       taxon=as.character(tax.info$TAXON), LM.a=tax.info$LM.a, LM.b=tax.info$LM.b, LM.p.ash=tax.info$LM.p.ash, p.b=tax.info$p.b,
                       temp.corr=temp.corr.sfpb, boot.num=boot.num)
          sink(paste("Output/", site, "_", gsub("/", "-", first.date), "_", gsub("/", "-", last.date), "_TALUS.CLIFF_", taxa.levels[a], ".txt", sep=""))
          print(out.tal)
          sink()
        }
      }
      else print(paste("Error: no method of production is specified for this taxon: ", taxa.levels[a]));Nintboots.full[[a]] = NA;Bintboots.full[[a]]= NA;Pintboots.full[[a]] = NA
      
      if (habitat == "COBBLE" & dim(tax.info)[1] != 0){
        Pboots.cob[,a] <- out$Pboots
        names(Pboots.cob)[a] <- as.character(tax.info$TAXON)
        Bboots.cob[,a] <- out$Bboots
        names(Bboots.cob)[a] <- as.character(tax.info$TAXON)
        Nboots.cob[,a] <- out$Nboots
        names(Nboots.cob)[a] <- as.character(tax.info$TAXON)
        Nintboots.out <- out$Nintboots
        Nintboots.cob[is.na(Nintboots.out) | is.null(Nintboots.out)]
        if (is.null(Nintboots.out) | all(is.na(Nintboots.out)) | all(is.nan(Nintboots.out))){
          Nintboots.cob[a,] <- 0
        }  
        else if(all(Nintboots.out == 0)){
          Nintboots.cob[a,] <- 0
        } 
        else {
          Nintboots.cob[a,] = apply(data.frame(Nintboots.out), 2, mean, na.rm = T)
        } 
        #Pintboots.out[which(Pintboots.out == 0)] = NA
        #browser()
        rownames(Nintboots.cob)[a] = as.character(tax.info$TAXON)
        Pintboots.out <- out$Pintboots
        Pintboots.out[is.na(Pintboots.out) | is.null(Pintboots.out)] <- 0
        if (is.null(Pintboots.out) |all(is.na(Pintboots.out)) | all(is.nan(Pintboots.out))){
          Pintboots.cob[a,] <- 0
        }  
        else if(all(Pintboots.out == 0)){
          Pintboots.cob[a,] <- 0
        } 
        else {
          Pintboots.cob[a,] = apply(data.frame(Pintboots.out), 2, mean, na.rm = T)
        } 
        #Pintboots.out[which(Pintboots.out == 0)] = NA
        rownames(Pintboots.cob)[a] = as.character(tax.info$TAXON)
        Bintboots.out <- out$Bdateboots
        if(is.null(Bintboots.out) | all(is.na(Bintboots.out))){
          Bintboots.cob[a,] <- 0
        }
        else if(all(Bintboots.out == 0)){
          Bintboots.cob[a,] <- 0
        }
        else {
          Bintboots.cob[a,] = apply(data.frame(Bintboots.out), 2, mean, na.rm = T)
        }
        #Bintboots.out[which(Bintboots.out == 0)] = NA
        #Bintboots.cob[a,] = apply(data.frame(Bintboots.out), 2, mean, na.rm = T)
        rownames(Bintboots.cob)[a] = as.character(tax.info$TAXON)
      
       
       if (is.null(Nintboots.out) | all(is.na(Nintboots.out)) | all(is.nan(Nintboots.out))){
         Nintboots.full[[a]] <- Bintboots.full[[a]] <- Pintboots.full[[a]] <- 0
       }  
       else if(all(Nintboots.out == 0)){
         Nintboots.full[[a]] <- Bintboots.full[[a]] <- Pintboots.full[[a]] <- 0
       } 
       else {
         # browser()
         Nintboots.full[[a]] = Nintboots.out;Bintboots.full[[a]]= Bintboots.out;Pintboots.full[[a]] = Pintboots.out;
         # lapply(list(Nintboots.full,Bintboots.full,Pintboots.full), FUN = function(x){ names(x[a]) <- as.character(tax.info$TAXON)
         # return(x)})
       }
      }
       
      else if (habitat == "DEPOSITIONAL" & dim(tax.info)[1] != 0){
        Pboots.dep[,a] <- out$Pboots
        names(Pboots.dep)[a] <- as.character(tax.info$TAXON)
        Bboots.dep[,a] <- out$Bboots
        names(Bboots.dep)[a] <- as.character(tax.info$TAXON)
        Nboots.dep[,a] <- out$Nboots
        names(Nboots.dep)[a] <- as.character(tax.info$TAXON)
      }
      else if (habitat == "TALUS.CLIFF" & dim(tax.info)[1] != 0){
        Pboots.tal[,a] <- out$Pboots
        names(Pboots.tal)[a] <- as.character(tax.info$TAXON)
        Bboots.tal[,a] <- out$Bboots
        names(Bboots.tal)[a] <- as.character(tax.info$TAXON)
        Nboots.tal[,a] <- out$Nboots
        names(Nboots.tal)[a] <- as.character(tax.info$TAXON)
      }
      else if (habitat == "ALL" & dim(tax.info)[1] != 0){
        Pboots.cob[,a] <- out.cob$Pboots
        Pboots.dep[,a] <- out.dep$Pboots
        Pboots.tal[,a] <- out.tal$Pboots
        Bboots.cob[,a] <- out.cob$Bboots
        Bboots.dep[,a] <- out.dep$Bboots
        Bboots.tal[,a] <- out.tal$Bboots
        Nboots.cob[,a] <- out.cob$Nboots
        Nboots.dep[,a] <- out.dep$Nboots
        Nboots.tal[,a] <- out.tal$Nboots
        names(Pboots.cob)[a] <- names(Pboots.dep)[a] <- names(Pboots.tal)[a] <- names(Bboots.cob)[a] <- names(Bboots.dep)[a] <- names(Bboots.tal)[a] <- names(Nboots.cob)[a] <- names(Nboots.dep)[a] <- names(Nboots.tal)[a] <- as.character(tax.info$TAXON)
      }
    }
    # browser()
        full_lists = list(Nintboots.full,Bintboots.full,Pintboots.full)
    full_lists = lapply(full_lists, function(x){names(x) <- as.character(taxa.levels)
      return(x)})
    list(Pboots.cob=Pboots.cob, Pboots.dep=Pboots.dep, Pboots.tal=Pboots.tal, Bboots.cob=Bboots.cob, Bboots.dep=Bboots.dep, Bboots.tal=Bboots.tal, Nboots.cob=Nboots.cob, Nboots.dep=Nboots.dep, Nboots.tal=Nboots.tal, Nintboots.cob = Nintboots.cob, Pintboots.cob = Pintboots.cob, Bintboots.cob = Bintboots.cob,
         Pintboots.full = full_lists[[3]], Bintboots.full = full_lists[[2]], Nintboots.full=full_lists[[1]])
    }
 }

