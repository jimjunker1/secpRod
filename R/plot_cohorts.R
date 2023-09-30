#' @name plot_cohorts
#' @title plot_cohorts
#' @description This function plots size-frequency data of a single taxon over time
#' @param taxaSampleListMass a data.frame of long format returned from \code{convert_length_to_mass()} function
#' @param param character. A string of 'length' or 'mass' that describes what measurement should be plotted
#' @param massClass character. The column name of the mass measurement (e.g., afdm_mg, dm_mg, etc.)
#' @param ... additional arguments passed to function
#' @return returns a histogram of the plot of the relative frequency of size or mass classes for a single taxon for all sampling dates
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom rlang .data
#' @importFrom graphics axis
#' @importFrom graphics grconvertX
#' @importFrom graphics hist
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics segments
#' @export

plot_cohorts = function(taxaSampleListMass = NULL, param = c('length','mass'), massClass = 'afdm_mg',...){
  # require(purrr)
  #### Tests ####
  if(length(unique(taxaSampleListMass$taxonID)) > 1) stop("Error: The number of taxa is >1. Only pass single taxon at a time.")
  if(any(taxaSampleListMass$n_m2 <0)) stop("Error: some density values < 0. Only positive values accepted.")
  # if(relFreq && max(n_m2) >1) stop("Error: `relFreq` is TRUE and the density values exceed 1. Check if data are relative frequencies or set `relFreq` to FALSE.")
  #### End Tests ####
  taxonID = unique(taxaSampleListMass$taxonID)

  # remove any length classes with all zeros
  allZeros = taxaSampleListMass %>%
    group_by(.data$lengthClass) %>%
    dplyr::summarise(n = sum(.data$n_m2, na.rm = TRUE)) %>%
    dplyr::filter(.data$n > 0) %>%
    dplyr::select(.data$lengthClass) %>%
    unlist

  taxaSampleListMass = taxaSampleListMass %>% dplyr::filter(.data$lengthClass %in% allZeros)

  ## Are the samples still separated by repID
    if(any(grepl("repID",names(taxaSampleListMass)))){
    repForm = as.formula(paste0("n_m2 ~ dateID + lengthClass + ",eval(massClass)))
      df = aggregate(repForm, data = taxaSampleListMass, FUN = 'sum')
    } else{
      df = taxaSampleListMass
    }
    # split by sampling dates and create length or mass vectors
  # if(param == 'length'){.
  #   plotDf = df[,c("dateID", "lengthClass", "n_m2")]
  # } else{
  #   plotDf = df[,c("dateID", eval(massClass),"n_m2")]
  # }
  plotDf = df[,c('dateID','lengthClass',eval(massClass), 'n_m2')]
  plotDf$n_m2 <- round(plotDf$n_m2)# this currently rounds to nearest integer. This is probably okay for data with high densities, but will need to change to account for low densities in next iteration.
  plotDfDateSplit = split(plotDf, as.factor(plotDf$dateID))
  plotDfDateSplitProps = lapply(plotDfDateSplit, function(x){
   x[rep(row.names(x), times = x$n_m2),]
  })
  plotDfExpanded = do.call('rbind', plotDfDateSplitProps)
  plotDfExpanded = plotDfExpanded[,c('dateID','lengthClass',eval(massClass))]

# set plotting parameters in base graphics
  # determine maximum x  y
  if(param == 'length'){
    xMax = max(df['lengthClass'])

    } else if(param == 'mass'){
      xMax = max(df[massClass])
    }
  # determine sampling dates
  nDates = length(unique(df$dateID))
  dateVec = sort(unlist(unique(df$dateID)))
  # determine maximum y
  parStore <- vector(mode = 'numeric', length = nDates)
  brks <- seq(0,(ceiling(xMax)+1), by = 1)
  for(j in 1:nDates){
    dateDf <- with(.data, subset(plotDfExpanded, dateID == dateVec[j]))
    if(param == 'length'){
      lVec <- unlist(dateDf['lengthClass'])
    } else{
      lVec <- unlist(dateDf[eval(massClass)])
    }
    if(length(lVec) == 0){
      parStore[j] <- NA_real_
    } else{
      # browser()
      h = hist(lVec, breaks = brks, plot = FALSE)
      h$counts <- h$counts/sum(h$counts)
      parStore[j] = max(h$counts)
      }
  }
  yMax = max(parStore, na.rm = TRUE)

  # set plot columns

  pfun <- function(){
  if(nDates < 6){
    par(mfcol = c(nDates,1))
  } else if(nDates%%2>0){
    par(mfcol = c((nDates+1)/2,2))
  } else{
    par(mfcol = c(nDates/2,2))
  }

  par(cex = 0.6,las=1,xaxs="i")
  par(mar = c(0, 1, 0, 0), oma = c(4, 4, 2, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))

  xlimcoord<-NA

  for(i in 1:nDates){
    dateDf <- with(.data, subset(plotDfExpanded, dateID == dateVec[i]))
    if(param == 'length'){
    lVec <- unlist(dateDf['lengthClass'])
    } else{
      lVec <- unlist(dateDf[eval(massClass)])
    }
    if(length(lVec) == 0){
      plot.new()
    } else{
      h = hist(lVec, breaks = brks, col = "lightgrey", border = "grey", main = "", xlim = c(0, xMax), probability = TRUE, axes=F,las=1,ylim=c(0,yMax), plot = FALSE)
      h$counts <- h$counts/sum(h$counts)
      plot(h, freq=TRUE,col = "lightgrey", border = "grey", main = "", xlim = c(0, xMax), axes=F,las=1,ylim=c(0,yMax))
    }

    xlimcoord[i]<-grconvertX(x=xMax, "user", "ndc")

    mtext(dateVec[i], side = 3, line = -1, adj = 0.1, cex = 0.6,col = "grey40")

    if(nDates<8){
      if (i %in% c((nDates/2),nDates))
        axis(1, col = "grey40", col.axis = "grey20", at = seq(0,round(xMax,3),length.out=10))

      axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(yMax,2),(round(yMax/3,2)) ))


    }else if(nDates%%2>0){
      if (i %in% c(((nDates+1)/2),nDates))
        axis(1, col = "grey40", col.axis = "grey20", at = seq(0,round(xMax,3),length.out =10))
      if (i %in% 1:((nDates+1)/2))
        axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(yMax,2),(round(yMax/3,2)) ))

    }else{

      if (i %in% c((nDates/2),nDates))
        axis(1, col = "grey40", col.axis = "grey20", at = seq(0,round(xMax,3),length.out=10))
      if (i %in% 1:(nDates/2))
        axis(2, col = "grey40", col.axis = "grey20", at = seq(0,round(yMax,2),(round(yMax/3,2)) ))
    }
  }
  if(param == 'length'){
  mtext("Length (mm)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
        col = "grey20")
  mtext("Proportion", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
        col = "grey20",las=0)
  } else{
    mtext("Mass (mg)", side = 1, outer = TRUE, cex = 0.7, line = 2.2,
          col = "grey20")
    mtext("Proportion", side = 2, outer = TRUE, cex = 0.7, line = 2.2,
          col = "grey20",las=0)
  }
  mtext(eval(taxonID), size = 3, outer = TRUE)
  }
  # debugonce(pfun)
  return(pfun())

  }
