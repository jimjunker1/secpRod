# ==============================================================================
# ==============================================================================
#     SCRIPT TO IMPLEMENT ALGORITHM FOR SOLVING FOR GROWTH PARAMETERS
#     VIA THE METHOD OUTLINED IN LLOYD-JONES ET AL. USES AN MM ALGORI
#     THM AND THE OPTIM FUNCTION TO OPTIMISE OVER LENGTH FREQUENCY DA
#     TA SETS.
# ==============================================================================
# ==============================================================================
# Remove any objects to clear the slate
rm(list = ls())
# ------------------------------------------------------------------------------
# Source the function files needed
# ------------------------------------------------------------------------------
# setwd("~/Desktop/bsc_r_codes_data/") # Set this path
here::i_am("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_main.R")
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_log_likelihood_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_mean_length_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_mean_var_optim_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_plot_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_seas_integral_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_seas_root_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_var_ricker_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_pi_calc_func.R"))
# ------------------------------------------------------------------------------
# Data preliminaries
# ------------------------------------------------------------------------------
# Read in the data set on the asymptotic males that was gathered through pots
lfd.big.males.females <- read.table(here::here("ignore/readings for len-freq function/Zhou et al. supp/data/bsc/bsc_lfd_pot.txt"), header = T)
# Pull out the lengths and the dates from these data files
lfd.big.males.females.dates   <- as.Date(lfd.big.males.females$Date, "%d/%m/%y")
lfd.big.males.females.lengths <- lfd.big.males.females$Carapacewidth
# Pull out the year and month information from these dates
lfd.big.males.females.dates.year   <- format(lfd.big.males.females.dates, '%Y')
lfd.big.males.females.dates.months <- format(lfd.big.males.females.dates, '%m')
# Read in the trawl data on males that contains juvenile recruitment and adults
lfd.trawl.males.females <- read.table(here::here("ignore/readings for len-freq function/Zhou et al. supp/data/bsc/bsc_lfd_trawl.txt"))
# Pull out the lengths and the dates from these data files
lfd.trawl.males.females.dates   <- as.Date(lfd.trawl.males.females $Date,
									"%d/%m/%Y")
lfd.trawl.males.females.lengths <- lfd.trawl.males.females$Carapace.width
# If both pot and trawl data are to be used set wit.pot = 1
with.pot = 1
# Concatenate the necessary elements from each file into a common
# dates and lengths array
if (with.pot == 1)
{
  lfd.dates    <- c(lfd.trawl.males.females.dates,   lfd.big.males.females.dates)
  lfd.lengths  <- c(lfd.trawl.males.females.lengths, lfd.big.males.females.lengths)
  combined.sex <- c(lfd.trawl.males.females$Sex,     lfd.big.males.females$Sex)
} else
{
  lfd.dates    <- c(lfd.trawl.males.females.dates)
  lfd.lengths  <- c(lfd.trawl.males.females.lengths)
  combined.sex <- c(lfd.trawl.males.females$Sex)
}
# ------------------------------------------------------------------------------
# Subset for males or females
# ------------------------------------------------------------------------------
# For a combined run set comb.on = 1. To run just males set males.on to 1.
# To run females set both males.on and comb.on to 0.
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
no.grps.par     <- 3
mu.init.par     <- c(0.5, 185, 50, 70)
theta.const     <- 0
thetas.par      <- c(-0.2, 0.2)
#thetas.par      <- c(0.68338122,   0.05821895)
var.init.par    <- c(19, 0.015)
yrs.old.par     <- c(0, 1, 2)
str.mnth.par    <- 1
num.months.seq  <- seq(1, num.months)
# ------------------------------------------------------------------------------
# Initialise the parameters of the model
# ------------------------------------------------------------------------------
num.inds <- length(months)                     # Number of individuals
no.grps  <- no.grps.par
pi.init  <- rep(1 / no.grps, num.months)
pis      <- matrix(rep(pi.init, each = no.grps), nrow = num.months,
                   ncol = no.grps)
mu.init  <- mu.init.par
if (theta.const == 1)
{
	theta.1    <- thetas.par[1]
	theta.2    <- thetas.par[2]
	max.contr  <- (1 / (2 * pi)) *
                  acos(theta.1  /
                  (sqrt(theta.2 ^ 2 +
               	  theta.1 ^ 2)))                # Calculates max of seas curve
	theta.2    <- (theta.1 * (sqrt(1 - cos(2 *
              	  pi * max.contr) ^ 2))) /
                  cos(2 * pi * max.contr)       # Theta 2 constrained by max
} else
{
	theta.1    <- thetas.par[1]
	theta.2    <- thetas.par[2]
}
k0      <- mu.init[1]
linf    <- mu.init[2]
mu.yr   <- mu.init[3:length(mu.init)]
thetas  <- c(theta.1, theta.2)
var.par <- var.init.par                       # Variance fun parameter vector
# ------------------------------------------------------------------------------
# Initialise the likelihood and set tolerence
# ------------------------------------------------------------------------------
mean.mnth.coh.mn <<- matrix(0, nrow = num.months, ncol = no.grps)
var.mnth.coh.mn  <<- matrix(0, nrow = num.months, ncol = no.grps)
for (i in seq(1, no.grps))
{
  mean.mnth.coh.mn[, i] <-  sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1,
                      	         theta.2 = theta.2, linf = linf , mu.yr = mu.yr,
                      	         yrs.old = yrs.old.par[i], str.mnth = str.mnth.par)
}
# Calculate the variances given the current update of the parameters
for (i in seq(1, no.grps))
{
  var.mnth.coh.mn[, i] <- sapply(mean.mnth.coh.mn[, i], BscVar, var.par.1 = var.par[1],
                                var.par.2 = var.par[2])
}
if (theta.const == 1)
{
	pars <- c(k0, linf, thetas[1], var.par, mu.yr)
} else {
	pars <- c(k0, linf, thetas, var.par, mu.yr)
}
MeanVarOptim(pars)
log.like.full <- -10e5
tol           <- 10e-6
log.like.old  <- -Inf
# ------------------------------------------------------------------------------
# Run while loop over procedure until convergence
# -----------------------------------------------
while (log.like.full - log.like.old > tol) {

  # Shift the current likelihood to the old likelihood

  log.like.old <- log.like.full


  # Calculate the pi for each group in each month
  # ---------------------------------------------

  # Returns a vector of pi with each column representing a month
  # and each row a group. Row 1 the largest. Row 2 the yr olds
  # and row 3 the juveniles

  pis <- t(sapply(num.months.seq, PiCalc))


  # Optimise the parameters for the means
  # -------------------------------------

  # Initialise and optimise

  if (theta.const == 1)
  {
	pars <- c(k0, linf, thetas[1], var.par, mu.yr)
  } else {
	pars <- c(k0, linf, thetas, var.par, mu.yr)
  }
  optim.means.var <- optim(pars, MeanVarOptim, control = list(maxit = 100000))
  pars            <- optim.means.var$par
  # Ask if optim converged

  print("Did optim converge?")
  print(optim.means.var$convergence)

  # Re-define the global parameters
  if (theta.const == 1)
  {
  	k0         <- pars[1]
  	linf       <- pars[2]
  	theta.1    <- pars[3]
  	var.par.1  <- pars[4]
  	var.par.2  <- pars[5]
  	mu.yr      <- pars[6:length(pars)]
  } else
  {
  	k0         <- pars[1]
  	linf       <- pars[2]
  	theta.1    <- pars[3]
  	theta.2    <- pars[4]
  	var.par.1  <- pars[5]
  	var.par.2  <- pars[6]
  	mu.yr      <- pars[7:length(pars)]
  }

  # If male or female we keep thetas fixed so turn off thetas
  # above and turn those on below. Look in bsc_mean_var_func.R
  # for more details
  if (theta.const == 1)
  {
  	theta.2   <- (theta.1 * (sqrt(1 - cos(2 * pi * max.contr)^2))) /
                 cos(2 * pi * max.contr)
  }

  # Calculate the means again for the final likelihood update
  for (i in seq(1, no.grps))
  {
  	mean.mnth.coh.mn[, i] <- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1,
                      	            theta.2 = theta.2, linf = linf , mu.yr = mu.yr,
                      	            yrs.old = yrs.old.par[i], str.mnth = str.mnth.par)
  }
  # Calculate the variances given the current update of the parameters
  for (i in seq(1, no.grps))
  {
  	var.mnth.coh.mn[, i] <- sapply(mean.mnth.coh.mn[, i], BscVar, var.par.1 = var.par.1,
                                   var.par.2 = var.par.2)
  }


  # Evaluate the likelihood

  log.like.full <- sum(sapply(num.months.seq, LogLikelihood2))


  # Give a plot of the current state of the model versus the data

  # BscPlotNew(c(theta.const, pars))

  # Pars including sigma^2 linf and theta 2

  pars.inc <- c(pars, theta.2, BscVar(var.par.1, var.par.2, linf))
  var.par  <- c(var.par.1, var.par.2)
  thetas   <- c(theta.1, theta.2)

  # Print out the loglikelihood, tolerance, and parameters

  print(c(log.like.full, log.like.full - log.like.old))
  print(pars.inc)

  # Tally up parameters for AIC and BIC - k, theta1, theta2, linf, mu0s, pis, var par, means

  pars.num <- length(pars)
  pi.num   <- length(as.numeric(pis)) - length(pis[, dim(pis)[2]])
  no.pars  <- pars.num + pi.num
  AIC      <- 2 * no.pars - 2 * log.like.full
  BIC      <- -2 * log.like.full + no.pars * log(num.inds)
  print(paste0("Information criteria - AIC ", round(AIC, 2), " BIC ",  round(BIC, 2)))
}
# ----------------------------------------------------------------------------------------
# Post run summary
# ----------------------------------------------------------------------------------------
BscPlotNew(c(theta.const, pars))
if (theta.const == 1)
{
  names(pars.inc) <- c("k0", "linf", "theta1", "a", "b",
                       "mu1", "mu2", "theta2", "sigma_linf")
} else
{
 names(pars.inc) <- c("k0", "linf", "theta1", "theta2", "a", "b", "mu1",
                      "mu2", "theta2", "sigma_linf")
}
library(knitr)
pars.inc.2 <- t(as.data.frame(pars.inc))
colnames(pars.inc.2) <- names(pars.inc)
kable(data.frame(pars.inc.2, AIC, BIC), row.names = F)










