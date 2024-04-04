## LEN_FREQ, plotting length frequency histograms of invert data, Author: Jim Junker
## Last Updated: Sept. 2016
##
##Purpose: This code automates the plotting of length-frequency histograms of specified site, habitat, and taxon for a specified date range.
##	This function is meant to be used with a matrix of invertebrate size-classes(mm) and abundance (m^-2). There are five preceeding columns
##	SITE, SAMPLE, DATE, HABITAT, TAXON.

##	Date should maintain the following form with exact headings on columns 1-5

###	SITE	SAMPLE	DATE		HABITAT	TAXON		0.25	0.5	1	1.5
#	HVER	2		08/02/2011	COBBLE	OSTRACOD	0	0	2	2
#	HVER	2		08/02/2011	COBBLE	MIDGE		0	0	10	7
#	HVER	2		08/02/2011	COBBLE	OSTRACOD	0	0	6	10
#	HVER	2		08/02/2011	COBBLE	TUBIFICID	0	10	5	0
#	HVER	4		08/02/2011	COBBLE	LIMNOPHORA	0	0	0	0
#	.
#	.
#	.

## Notes:
##	SITE can handle an unlimited number of sites
##	SAMPLE number is detected and can handle unequal sample numbers across site/date/habitat/etc.
##	DATE This is the date of sampling. This must be in the form of MM/DD/YYY (sample date: "1/5/2008", "9/15/2007", etc.). This maintains the format for further
##		secondary production estimation.
##	HABITAT For habitat specific production. This is necessary for further secondary production estimation. Can have unlimited number of habitat types
##	TAXON This is the taxon or cohort information used for each group to summarize. This code will summarize the taxon with exactly matching names within samples.
##	Column 6 - n (as many as needed) These headings denote mm size classes. The values in these columns are abundances sample (m^2).
##		There is no requirement for these to be in order or in constant intervals. All columns should be filled in with counts or zeros. Truly missing data
##		should be filled with "NA". There is no need to include all non-occurrences of taxa, but will work fine if it does.
##
##	When loading in the data file be sure to turn name checking off otherwise
##	(i.e. DATA <- read.table(C:\DATA.txt, header = T, sep = "\t", quote = "", strip.white = T, row.names = F, check.names = F))
##

##Required packages 'plyr', 'dplyr', 'tidyr', 'reshape2', 'ggplot2', 'chron', 'grid', 'gridExtra'
if(!require("pacman")) install.packages("pacman")
library(pacman)
package.list <- c("plyr", "dplyr", "tidyr", "reshape2", "ggplot2", "chron", "grid", "gridExtra")
p_load(char = package.list, install = T)

rm("package.list")

####

len_freq <- function(DATA, site, TAXA, habitat, first.date, last.date, fun, ...){
theme_set(theme_classic())

graphics.off()

  DATA_SUM <- aggregate(DATA[c(names(DATA[7:(dim(DATA)[2])]))], by = DATA[c("SITE", "DATE", "JULIAN", "HABITAT", "TAXON")], FUN = fun)
  DATA_COL <- DATA_SUM[,1:5]
  DATA_PROP <- t(apply(DATA_SUM[6:(dim(DATA_SUM)[2])],1, prop.table))
  PROP_TABLE <- data.frame(DATA_COL, DATA_PROP, check.names = F)

#source("C:/Users/Jim/Documents/Projects/Iceland/Bug Samples/Secondary Production/Secondary Production R code suite/len_freq/len_freq_plot.txt")
if(missing(site) & missing(TAXA) & missing(habitat) & missing(first.date) & missing(last.date)) {

  sites = levels(DATA$SITE)
  hab = levels(DATA$HABITAT)
  taxon = levels(DATA$TAXON)

   for(i in sites){
     data_site = PROP_TABLE[which(PROP_TABLE$SITE == i),]
   for(j in hab){
    data_hab = data_site[which(data_site$HABITAT == j),]
   for(k in taxon){

  data_tax = data_hab[which(data_hab$TAXON == k),]
	dates = as.character(sort(unique(data_tax$JULIAN)))
	pltList = list()

   for(l in dates){
	data_loop = data_tax[which(data_tax$JULIAN == l),];
	date_data = gather(data_loop, size, rel.freq, 6:(dim(data_loop)[2]));
	v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
	ymax = max(date_data$rel.freq) + 0.05
	xmax = max(date_data$size[which(date_data$rel.freq != 0)])

	pltList[[l]] <- ggplot(date_data, aes( x = as.numeric(size), y = rel.freq)) +
	geom_bar(stat = "identity", width = 0.99) +
	scale_y_continuous(limits = c(0, as.numeric(ymax)), expand = c(0,0)) +
	scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
	geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
	theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == l)])))

  args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = "Relative Frequency", bottom = "Size (mm)")
  pdf(file = paste(data_tax$SITE,"_",data_tax$TAXON,"_plot.pdf"))
  do.call(grid.arrange, args.list)
  dev.off()
    }
   }
  }
 }
} #else if(!missing(site) & missing(TAXA) & missing(habitat) & missing(first.date) & missing(last.date)) {
  else if(!missing(site)) {
 data_site = PROP_TABLE[which(PROP_TABLE$SITE == site),]

 hab = levels(data_site$HABITAT)
 taxon = levels(data_site$TAXON)

  for(j in hab){

   data_hab = data_site[which(data_site$HABITAT == habitat)]

    for(k in taxon){
      data_tax = data_hab[which(data_hab$TAXON == k),]
      dates = as.character(sort(unique(data_tax$JULIAN)))
      pltList = list()

      for(l in dates){
        data_loop = data_tax[which(data_tax$JULIAN == l),];
        date_data = gather(data_loop, size, rel.freq, 6:(dim(data_loop)[2]));
        v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
        ymax = max(date_data$rel.freq) + 0.05
        xmax = max(date_data$size[which(date_data$rel.freq != 0)])

        pltList[[l]] <- ggplot(date_data, aes( x = as.numeric(size), y = rel.freq)) +
          geom_bar(stat = "identity", width = 0.99) +
          scale_y_continuous(limits = c(0, as.numeric(ymax)), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
          geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
          ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == l)])))

          args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = "Relative Frequency", bottom = "Size (mm)")
          pdf(file = paste(data_tax$SITE,"_",data_tax$TAXON,"_plot.pdf"))
          do.call(grid.arrange, args.list)
          dev.off()
      }
    }
  }
 } #else if(!missing(site) & !missing(habitat) & missing(TAXA) & missing(first.date) & missing(last.date)) {
else if(!missing(site) & !missing(habitat)) {

data_site_hab = PROP_TABLE[which(PROP_TABLE$SITE == site) & which(PROP_TABLE$HABITAT == habitat),]
taxon = levels(data_site_hab$TAXON)

  for(k in taxon){
    data_tax = data_site_hab[which(PROP_TABLE$TAXON == k),]
    dates = as.character(sort(unique(data_tax$JULIAN)))
    pltList = list()

    for(l in dates){
      data_loop = data_tax[which(data_tax$JULIAN == l),];
      date_data = gather(data_loop, size, rel.freq, 6:(dim(data_loop)[2]));
      v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
      ymax = max(date_data$rel.freq) + 0.05
      xmax = max(date_data$size[which(date_data$rel.freq != 0)])

      pltList[[l]] <- ggplot(date_data, aes( x = as.numeric(size), y = rel.freq)) +
        geom_bar(stat = "identity", width = 0.99) +
        scale_y_continuous(limits = c(0, as.numeric(ymax)), expand = c(0,0)) +
        scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
        geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == l)])))

        args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = "Relative Frequency", bottom = "Size (mm)")
        pdf(file = paste(data_tax$SITE,"_",data_tax$TAXON,"_plot.pdf"))
        do.call(grid.arrange, args.list)
        dev.off()
    }
  }
} #else if(!missing(site) & !missing(habitat) & !missing(TAXA) & missing(first.date) & missing(last.date)){
else if(!missing(site) & !missing(habitat) & !missing(TAXA)){

  data_site_hab_tax = PROP_TABLE[which(PROP_TABLE$SITE == site) & which(PROP_TABLE$HABITAT == habtita) & which(PROP_TABLE$TAXON == TAXA),]
  dates = as.character(sort(unique(data_sub$JULIAN)))
      pltList = list()

      for(l in dates){
        data_loop = data_site_hab_tax[which(data_site_hab_tax$JULIAN == l),];
        date_data = gather(data_loop, size, rel.freq, 6:(dim(data_loop)[2]));
        v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
        ymax = max(date_data$rel.freq) + 0.05
        xmax = max(date_data$size[which(date_data$rel.freq != 0)])

        pltList[[l]] <- ggplot(date_data, aes( x = as.numeric(size), y = rel.freq)) +
          geom_bar(stat = "identity", width = 0.99) +
          scale_y_continuous(limits = c(0, as.numeric(ymax)), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
          geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
          ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == l)])))

          args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = "Relative Frequency", bottom = "Size (mm)")
          pdf(file = paste(data_tax$SITE,"_",data_tax$TAXON,"_plot.pdf"))
          do.call(grid.arrange, args.list)
          dev.off()
   }
  }
}
