# ==============================================================================
# Author: Luke Lloyd-Jones
# Date started: 15/12/2015
# Date updated: 29/12/2015    
# Script to do regular investigation of the initial data                                                  
# ==============================================================================
setwd("~/Dropbox/AAUni/Blue_Swimmer/biometrics_review/r_codes_data/")
#install.packages("mixtools")
install.packages("mclust")
library(mclust)
# ------------------------------------------------------------------------------
# Data preliminaries
# ------------------------------------------------------------------------------
setwd("~/Dropbox/AAUni/Blue_Swimmer/biometrics_review/r_codes_data/")
# ------------------------------------------------------------------------------
# POTS
# ------------------------------------------------------------------------------
# Read in the data set on the asymptotic males that was gathered through pots
lfd.big.males.females <- read.table("Data/BSC/bsc_lfd_pot.txt", header = T)
# Pull out the lengths and the dates from these data files
lfd.big.males.females.dates   <- as.Date(lfd.big.males.females$Date, "%d/%m/%y")
lfd.big.males.females.lengths <- lfd.big.males.females$Carapacewidth
# Pull out the year and month information from these dates
lfd.big.males.females.dates.year   <- format(lfd.big.males.females.dates, '%Y')
lfd.big.males.females.dates.months <- format(lfd.big.males.females.dates, '%m')
# ------------------------------------------------------------------------------
# TRAWL
# ------------------------------------------------------------------------------
# Read in the trawl data on males that contains juvenile recruitment and adults
lfd.trawl.males.females <- read.table("Data/BSC/bsc_lfd_trawl.txt")
# Pull out the lengths and the dates from these data files
lfd.trawl.males.females.dates   <- as.Date(lfd.trawl.males.females $Date, 
									"%d/%m/%Y")
lfd.trawl.males.females.lengths <- lfd.trawl.males.females$Carapace.width
with.pot = 1
just.pot = 0
# Concatenate the necessary elements from each file into a common 
# dates and lengths array
if (with.pot == 1)
{
  lfd.dates   <- c(lfd.trawl.males.females.dates,   lfd.big.males.females.dates)
  lfd.lengths <- c(lfd.trawl.males.females.lengths, lfd.big.males.females.lengths)
  combined.sex <- c(lfd.trawl.males.females$Sex,    lfd.big.males.females$Sex)
} else if (just.pot == 1)
{
  lfd.dates    <- c(lfd.big.males.females.dates)
  lfd.lengths  <- c(lfd.big.males.females.lengths)
  combined.sex <- c(lfd.big.males.females$Sex)
} else
{
  lfd.dates   <- c(lfd.trawl.males.females.dates)
  lfd.lengths <- c(lfd.trawl.males.females.lengths)	
  combined.sex <- c(lfd.trawl.males.females$Sex)
}
# ------------------------------------------------------------------------------
# Subset for males or females
# ------------------------------------------------------------------------------
males.on <- 0
comb.on  <- 1
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
str.month.yr.1 <- 1
end.month.yr.1 <- 12
str.month.yr.2 <- 1
end.month.yr.2 <- 12
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
num.mnts        <- length(names(table(months)))
months.lst      <- as.numeric(names(table(months)))
# ------------------------------------------------------------------------------
# Fit regular mixtures to the data
# ------------------------------------------------------------------------------
library(ggplot2)
library(grid)
# Declare some plotting parameters
x.vals <- 0 : 200
x.lab.name <- c("Jan., 1985", "Feb., 1985", "Mar., 1985", "Apr., 1985", "May, 1985",
                "Jun., 1985", "Jul., 1985", "Aug., 1985", "Sep., 1985", "Oct., 1985",
                "Nov., 1985", "Dec., 1985", "Jan., 1986",
                "Feb., 1986", "Mar., 1986", "Apr., 1986", "May, 1986", "Jun., 1986")
# Set up the plotting parameters
dev.new(width = 15, height = 17/1.618, units = "cm")
plot.nrow <- 3
plot.ncol <- 6
pushViewport(viewport(layout = grid.layout(plot.nrow, plot.ncol)))
row.seq <- rep(seq(1, plot.nrow), each  = plot.ncol)
col.seq <- rep(seq(1, plot.ncol), times = plot.nrow)
n.clst <- 3
# Store the means and variances for later plotting
pi.clst    <- matrix(nrow = num.mnts, ncol = n.clst)
means.clst <- matrix(nrow = num.mnts, ncol = n.clst)
var.clst   <- matrix(nrow = num.mnts, ncol = n.clst)
bic        <- matrix(nrow = num.mnts, ncol = 1)
# Plot of data and mixture fits
for (mm in 1 : num.mnts) {
  month.lengths <- lengths[which(months == months.lst[mm])]
  lengths.sub.df <- data.frame(month.lengths)
  m <- ggplot(lengths.sub.df, aes(x = month.lengths)) +
       geom_histogram(aes(y = ..density..), 
       fill = rgb(0.5, 0.5, 0.5, 0.5),     
       binwidth = 7) + ylim(c(0, 0.05)) +
       xlab(x.lab.name[mm]) 
  #lengths.mm.clst <- Mclust(month.lengths, G = 3)
  lengths.mm.clst <- Mclust(month.lengths, G = n.clst)
  print(summary(lengths.mm.clst, parameters = TRUE))
  #plot(lengths.mm.clst)
  pis <- summary(lengths.mm.clst, parameters = TRUE)$pro
  means <- summary(lengths.mm.clst, parameters = TRUE)$mean
  vars  <- summary(lengths.mm.clst, parameters = TRUE)$var
  pi.clst[mm, ]    <- pis
  means.clst[mm, ] <- means
  var.clst[mm, ]   <- vars
  bic[mm, ]        <- summary(lengths.mm.clst, parameters = TRUE)$bic
  den.mu2 <- pis[1] * dnorm(x.vals, means[1], sqrt(vars[1]))
  den.mu1 <- pis[2] * dnorm(x.vals, means[2], sqrt(vars[2]))
  den.mu0 <- pis[3] * dnorm(x.vals, means[3], sqrt(vars[3]))
  x.df.mu2 <- data.frame(x.vals, den.mu2)
  x.df.mu1 <- data.frame(x.vals, den.mu1)
  x.df.mu0 <- data.frame(x.vals, den.mu0)
  if (n.clst == 3)
  {
    total.mixture <- den.mu2 + den.mu1 + den.mu0
  } else {
  	total.mixture <- den.mu2 + den.mu1 
  }
  x.vals.tot    <- c(x.vals, x.vals, x.vals)
  tot.df        <- data.frame(x.vals.tot, total.mixture)
  plot.mm  <- m +  
          geom_line(data = x.df.mu2, aes(x = x.vals, y = den.mu2), 
          linetype = 1, size = 1) +
          geom_line(data = x.df.mu1, aes(x = x.vals, y = den.mu1), 
          linetype = 2, size = 1) +
          geom_line(data = x.df.mu0, aes(x = x.vals, y = den.mu0), 
          linetype = 3, size = 1) +
          geom_line(data = tot.df,   aes(x = x.vals.tot, y = total.mixture), 
          linetype = 1, size = 0.5) +
          xlab(x.lab.name[mm]) 
  print(plot.mm, vp = viewport(layout.pos.row = row.seq[mm], 
        layout.pos.col = col.seq[mm]))
  #print(m, vp = viewport(layout.pos.row = row.seq[mm], 
  #      layout.pos.col = col.seq[mm]))
}
# ------------------------------------------------------------------------------
# Generate plots of the variances as a function of the mean
# With some fitted line through it
# ------------------------------------------------------------------------------
library(reshape2)
vc <- as.numeric(var.clst)
mc <- as.numeric(means.clst)
mn.vr <- data.frame(mc, vc)
c <- ggplot(mn.vr, aes(mc, vc))
c + stat_smooth(aes(size = 2)) + geom_point(aes(size = 2)) +
  xlab("Estimated mean (mm)") + ylab("Estimated variance") +
  theme(legend.position = "none")
# ------------------------------------------------------------------------------
# Generate plots of the month vs the mean for that month as 
# with some fitted line through the points
# ------------------------------------------------------------------------------
cbind(c(seq(1, 11), seq(1, 7)), means.clst)
mn.mnt.1 <- data.frame(seq(1, 33), as.numeric(means.clst[1:11, ]))
colnames(mn.mnt.1) <- c("Month", "Mean_len")
mn.mnt.2 <- data.frame(seq(1, 21), as.numeric(means.clst[12:18, ]))
colnames(mn.mnt.2) <- c("Month", "Mean_len")
mn.mnt <- rbind(mn.mnt.1, mn.mnt.2)
colnames(mn.mnt) <- c("Month", "Mean_len")
c <- ggplot(mn.mnt, aes(Month, Mean_len))
c + stat_smooth(aes(size = 2)) + geom_point(aes(size = 2)) +
  xlab("Month") + ylab("Estimated mean length (mm)") +
  theme(legend.position = "none")
# Use optim to solve for VB equation
k <- 0.5
linf <- 190
t0 <- 0.5
theta.init <- c(k, linf)
ls <- function(theta)
{
	k    <- theta[1]
	linf <- theta[2]
	#t0   <- theta[3]
	pd.lng <- linf * (1 - exp(-k * (mn.mnt[, 1] / 12 - t0)))
	# Sum of squares
	sum((mn.mnt[, 2] - pd.lng) ^2)
}
res <- optim(theta.init, ls)
theta  <- res$par
k <- theta[1]
linf <- theta[2]
t0   <- theta[3]
time <- seq(0, 3, 0.1)
pd.lng <- linf * (1 - exp(-k * (time - t0)))

# ------------------------------------------------------------------------------
# Table up the pis, means and the variances and put in 
# web resources
# ------------------------------------------------------------------------------
library(knitr)
clst <- round(cbind(pi.clst, means.clst, var.clst), 3)
colnames(clst) <- c("pi1", "pi2", "pi3", 
                    "mu1", "mu2", "mu3", 
                    "sigma_1", "sigma2", "sigma3")
rownames(clst) <- x.lab.name
kable(clst, format = "latex")
