# Mean function for calcuating the mean length for the current month
# ------------------------------------------------------------------


MeanLengthGomp <- function(month, k0, theta.1, theta.2 , linf, mu.yr, 
			               yrs.old, str.mnth) {
	
  # Computes the mean length of the distribution for the current month
  # based on an integration over the seasonal curve from a von Bertal
  # anffy growth model.
  #
  # Args:
  #  month:    Current month to calulate mean at. Numbered from Jan=0
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  linf:     Asymptotic length
  #  mu.yr:    Vector of starting means for each year modelled
  #  yrs.old:  How old the individuals in this group are
  #  str.mnth: The month that we start calculating means for 
  #
  # Returns:
  #  Mean length for the month in the current iteration
  
  # Test values 
  
  #month    = 0 
  #k0       = pars[1]
  #theta.1  = pars[2]
  #theta.2  = (theta.1 * (sqrt(1 - cos(2 * pi * max.contr)^2))) /
  #           cos(2 * pi * max.contr)  
  #linf     = pars[3]
  #mu.yr.1  = pars[4]
  #mu.yr.2  = pars[5]
  #yrs.old  = 0
  #str.mnth = 0
			  
  # Set each of the months to be the middle of the month (1/24 to be used
  # with the seasonal function.
  
  mm.val  <- (month %% 12) / 12 + 1 / 24 + yrs.old
  str.mid <- str.mnth / 12 + 1 / 24 + yrs.old
  end.mid <- 1 + 1 / 24 + yrs.old


  # Assess whether the parameters at this update cross the y=0 axis
  # This will allow us to assess whether we need to calculate roots 
  # for the seasonal function or not
  
  time.var  <- seq(0, 1, 0.01)
  seas.func <- SeasFunc(k0, theta.1, theta.2, time.var)
  is.neg.1  <- min(seas.func)
  is.neg.2  <- max(seas.func)
	  
  if (is.neg.1 < 0 & is.neg.2 > 0) {


	# Calculate the roots of the seasonal function 
	
    root.1 <- SeasRootCalc(k0, theta.1, theta.2, yrs.old)[1]
    root.2 <- SeasRootCalc(k0, theta.1, theta.2, yrs.old)[2]

	    
	# Given the roots calculate the integral over the previous years
	# This is equal to the integral up to the first root + the in
	
	int.root.1   <- SeasIntegFunc(k0, theta.1, theta.2, str.mid, root.1)
	int.root.2   <- SeasIntegFunc(k0, theta.1, theta.2, root.2,  end.mid)
	int.yrs.prvs <- yrs.old * (int.root.1 + int.root.2)
	
	    
    # Integral for those months less than root 1
    
    if (mm.val <= root.1) {
    	
      integral = int.yrs.prvs + 
                 SeasIntegFunc(k0, theta.1, theta.2, str.mid, mm.val)
    }
	    	
    # Integral for those months greater than root 1 but less than root 2

    if (mm.val > root.1 & mm.val < root.2) {
    	
    	integral = int.yrs.prvs +
                   SeasIntegFunc(k0, theta.1, theta.2, str.mid, root.1)
    }
	    
    # Integral for those months less than root 1
    
    if (mm.val >= root.2) {
    	
       int.root.1.1 <- SeasIntegFunc(k0, theta.1, theta.2, str.mid, root.1)
       int.root.2.1 <- SeasIntegFunc(k0, theta.1, theta.2, root.2,  mm.val)
	   integral     <- int.root.1.1 + int.root.2.1 + int.yrs.prvs
	  }
	    
    } else {
    
    # If the integral doesn't have a negative component
    int.tot   <- SeasIntegFunc(k0, theta.1, theta.2, 1/24, 1 + 1/24)
	int.yrs.prvs <- yrs.old * int.tot
    int.no.root <- SeasIntegFunc(k0, theta.1, theta.2, str.mid, mm.val)
    integral    <- int.no.root + int.yrs.prvs
    
    }
 
  #print(integral)
  
  # Depending on the year modelled the starting mean will be different
  # so we change the value returned depending what month we are calculating
  # the mean for
  # Create an index
  if (month == 0)
  {
  	mu.yr.ind <- 1
  } else
  {
    mu.yr.ind <- ceiling(month / 12)	 
  }
  return(exp(log(mu.yr[mu.yr.ind]) + (log(linf) - log(mu.yr[mu.yr.ind])) * (1 - exp(-integral))))
  #return(mu.yr[mu.yr.ind] + (linf - mu.yr[mu.yr.ind]) * (1 - exp(-integral)))
  
}
