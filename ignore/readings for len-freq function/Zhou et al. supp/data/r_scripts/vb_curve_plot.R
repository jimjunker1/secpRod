# ==============================================================================
# Author: Luke Lloyd-Jones
# Date started: 07/01/2016
# Date updated: 07/01/2016    
# Script to plot the von Bertalanffy growth curves for BSC results                                                  
# ==============================================================================
library(ggplot2)
library("reshape2")
# Estimates from model k, theta_1, theta_2, mu_inf
est.comb <- c(0.715, 0.315, 0.016,  179.4)
est.male <- c(0.565, 0.308, -0.473, 196.2)
est.fem  <- c(0.787, 0.0574, 0.435, 161.6)
# Set up a grid of time maybe three years
t <- seq(0.05, 6, 0.05)
# Length function
gvbgm <- function(theta, t)
{
	l0 <- 50
    k0 <- theta[1]
    str.time <- 0.05
    end.time <- t
    theta.1  <- theta[2]
    theta.2  <- theta[3]
    mu_inf   <- theta[4]
	zt <- k0 * (end.time - str.time) + 
           (theta.1 / (2 * pi)) * (sin(2 * pi * end.time) - sin(2 * pi * str.time)) - 
           (theta.2 / (2 * pi)) * (cos(2 * pi * end.time) - cos(2 * pi * str.time))
    lt <- l0 + (mu_inf - l0) * (1 - exp(-zt))
    return(lt)
}
# For each of the parameter sets calculate the mean curve over the time frame
# Combined
Combined <- gvbgm(est.comb, t)
# Male
Male <- gvbgm(est.male, t)
# Female
Female  <- gvbgm(est.fem, t)
# Make a data frame for interfacing with ggplot
mean.crvs <- data.frame(t, Combined, Male, Female)
# Melt it
mean.crvs.mlt <- melt(mean.crvs, id = "t")
# Plot it
ggplot(data = mean.crvs.mlt,
       aes(x = t, y = value, colour = variable)) +
       geom_line(size = 1) + ylim(c(40, 200)) + ylab("Length (mm)") + 
       labs(colour = "Analysis") + 
       xlab("Time (years)")
