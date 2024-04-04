# ==============================================================================
# Author: Luke Lloyd-Jones/Hien Nguyen
# Date started: 15/12/2015
# Date updated: 29/12/2015    
# Script to do goodness of fit test so that we can say it fits goood!                                                  
# ==============================================================================
# Remove any objects to clear the slate
rm(list = ls())
# ------------------------------------------------------------------------------
# Source the function files needed
# ------------------------------------------------------------------------------
setwd("~/Dropbox/AAUni/Blue_Swimmer/biometrics_review/r_codes_data/")
source("MTDGM_BSC/mtdgm_log_likelihood_func.R")
source("MTDGM_BSC/mtdgm_mean_length_func.R")
source("MTDGM_BSC/mtdgm_mean_var_optim_func.R")
source("MTDGM_BSC/mtdgm_plot_func.R")
source("MTDGM_BSC/mtdgm_seas_integral_func.R")
source("MTDGM_BSC/mtdgm_seas_root_func.R")
source("MTDGM_BSC/mtdgm_var_ricker_func.R")
#source("MTDGM_BSC/mtdgm_var_quad.R")
source("MTDGM_BSC/mtdgm_pi_calc_func.R")
# ------------------------------------------------------------------------------
# Data preliminaries
# ------------------------------------------------------------------------------
# Read in the data set on the asymptotic males that was gathered through pots
lfd.big.males.females <- read.table("Data/BSC/bsc_lfd_pot.txt", header = T)
# Pull out the lengths and the dates from these data files
lfd.big.males.females.dates   <- as.Date(lfd.big.males.females$Date, "%d/%m/%y")
lfd.big.males.females.lengths <- lfd.big.males.females$Carapacewidth
# Pull out the year and month information from these dates
lfd.big.males.females.dates.year   <- format(lfd.big.males.females.dates, '%Y')
lfd.big.males.females.dates.months <- format(lfd.big.males.females.dates, '%m')
# Read in the trawl data on males that contains juvenile recruitment and adults
lfd.trawl.males.females <- read.table("Data/BSC/bsc_lfd_trawl.txt")
# Pull out the lengths and the dates from these data files
lfd.trawl.males.females.dates   <- as.Date(lfd.trawl.males.females $Date, 
									"%d/%m/%Y")
lfd.trawl.males.females.lengths <- lfd.trawl.males.females$Carapace.width
with.pot = 1
# Concatenate the necessary elements from each file into a common 
# dates and lengths array
if (with.pot == 1)
{
  lfd.dates    <- c(lfd.trawl.males.females.dates,   lfd.big.males.females.dates)
  lfd.lengths  <- c(lfd.trawl.males.females.lengths, lfd.big.males.females.lengths)
  combined.sex <- c(lfd.trawl.males.females$Sex,    lfd.big.males.females$Sex)
} else
{
  lfd.dates    <- c(lfd.trawl.males.females.dates)
  lfd.lengths  <- c(lfd.trawl.males.females.lengths)	
  combined.sex <- c(lfd.trawl.males.females$Sex)
}
# ------------------------------------------------------------------------------
# Subset for males or females
# ------------------------------------------------------------------------------
males.on <- 0
comb.on  <- 0
if (males.on == 1 & comb.on != 1)
{
  # Males
  # -----
  males        <- which(combined.sex == 1)
  lfd.dates    <- lfd.dates[males]
  lfd.lengths  <- lfd.lengths[males]
} else if (males.on != 1 & comb.on != 1)
{
  # Females
  # -------
  females      <- which(combined.sex == 2)
  lfd.dates    <- lfd.dates[females]
  lfd.lengths  <- lfd.lengths[females]
}
# Pull out the year and month information from these dates
lfd.year   <- format(lfd.dates, '%Y')
lfd.months <- format(lfd.dates, '%m')
# Pull out the years and months that we are interested in 
# i.e., those that don't contain recruitment
str.month.yr.1 <- 2
end.month.yr.1 <- 8
str.month.yr.2 <- 2
end.month.yr.2 <- 5
lfd.85 <- which((lfd.year == '1985') & 
				  (as.numeric(lfd.months) %in% (str.month.yr.1:end.month.yr.1)))
lfd.86 <- which((lfd.year == '1986') &
				  (as.numeric(lfd.months) %in% (str.month.yr.2:end.month.yr.2)))
# ------------------------------------------------------------------------------
# Initialise the data for the model 
# ------------------------------------------------------------------------------	
lengths         <- lfd.lengths[c(lfd.85, lfd.86)]		
months.85       <- as.numeric(lfd.months[lfd.85]) - 1 # Jan = 0th month
months.86       <- as.numeric(lfd.months[lfd.86]) + 11
months          <- c(months.85, months.86)
num.months      <- length(names(table(months)))	
months.lst      <- as.numeric(names(table(months)))
# ------------------------------------------------------------------------------
# BIND up the data
# ------------------------------------------------------------------------------
comb.dat <- cbind(lengths, months)
dim(comb.dat) # 15065
save(comb.dat, file = "~/Desktop/comb_dat.Rdata")
mal.dat <- cbind(lengths, months)
dim(mal.dat) # 8845
save(mal.dat, file = "~/Desktop/mal_dat.Rdata")
fem.dat <- cbind(lengths, months)
dim(fem.dat) # 6220
save(fem.dat, file = "~/Desktop/fem_dat.Rdata")
# ------------------------------------------------------------------------------
# PIS, MEAN, AND VARIANCES
# ------------------------------------------------------------------------------
comb.list <- list(pis, mean.mnth.coh.mn, var.mnth.coh.mn)
names(comb.list) <- c("comb.pi", "comb.mean", "comb.var")
save(comb.list, file = "~/Desktop/comb_res_list.Rdata")

mal.list <- list(pis, mean.mnth.coh.mn, var.mnth.coh.mn)
names(mal.list) <- c("mal.pi", "mal.mean", "mal.var")
save(mal.list, file = "~/Desktop/male_res_list.Rdata")

fem.list <- list(pis, mean.mnth.coh.mn, var.mnth.coh.mn)
names(fem.list) <- c("fem.pi", "fem.mean", "fem.var")
save(fem.list, file = "~/Desktop/fem_res")
# ------------------------------------------------------------------------------
# Load up the lists
# ------------------------------------------------------------------------------
load("~/Desktop/comb_dat.Rdata")
load("~/Desktop/mal_dat.Rdata")
load("~/Desktop/fem_dat.Rdata")
load("~/Desktop/comb_res_list.Rdata")
load("~/Desktop/mal_res_list.Rdata")
load("~/Desktop/fem_res_list.Rdata")