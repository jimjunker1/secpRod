get_temp_dates <- function(DATA, site, habitat, first.date=NA, last.date=NA){
  if (habitat == "COBBLE" | habitat == "DEPOSITIONAL" | habitat =="TALUS.CLIFF"){
    if (is.na(first.date) | is.na(last.date)){
      data4 <- DATA[DATA$SITE == site & DATA$HABITAT == habitat,]
    }
    else{
      data4 <- DATA[DATA$SITE == site & DATA$JULIAN >= DATA$JULIAN[DATA$DATE == first.date][1] & DATA$JULIAN <= DATA$JULIAN[DATA$DATE == last.date][1] & DATA$HABITAT == habitat,]
    }
    data4 <- data4[!is.na(apply(data4[,7:(dim(data4)[2])], 1, sum)),]			#Remove rows with NA's (missing data)
    row.names(data4) <- NULL								#Reset row names to be sequential
    samp.julians <- data4 %>% group_by(DATE) %>%
      summarise(JULIAN = unique(JULIAN))
    # samp.julians <- sort(as.numeric(levels(as.factor(data4$JULIAN))))			#Define the sampling Julian dates
    # samp.dates <- numeric(length(samp.julians))
    # for (i in 1:length(samp.julians)){
    #   samp.dates[i] <- as.character(data4$DATE[data4$JULIAN == samp.julians[i]][1])	#Define the sampling dates
    # }
    # data.frame(DATE=samp.dates, JULIAN=samp.julians)
    data.frame(samp.julians)
  }
  else print("Error: invalid habitat")
}