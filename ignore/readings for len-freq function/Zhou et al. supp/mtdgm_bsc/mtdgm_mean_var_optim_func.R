# Function to optimise means using optim
# --------------------------------------
	
MeanVarOptim <- function(pars) {
  # Calculates the log likelihood 
  #
  # Args:
  #  pars: Vector containing the parameters to optimise
  #    k0:        Mean  k paramter for VB model
  #    thetas
  #    	theta.1:   Seasonality parameter 1
  #    	theta.2:   Seasonality parameter 1
  #    linf:      Asymptotic length
  #    mu.yr  :   vector containing all the strating months means
  #    var.pars.x: The variance parameters
  #    theta.const: boolean parametr that specifies if 1 that we should constrain
  #                 the variance function
  # Returns:
  #  The negative of the log likelihood for use in the OPTIM function
  
  # Declare each of the parameters to names unique to inside the fucntion
  # print(theta.const)
  if (theta.const == 1)
  {
  	k0.fun        <- pars[1]
  	linf.fun      <- pars[2]
  	theta.1.fun   <- pars[3]
  	var.par.1.fun <- pars[4]
  	var.par.2.fun <- pars[5]
  	mu.yr.fun     <- pars[6:length(pars)]
  } else {
  	k0.fun        <- pars[1]
  	linf.fun      <- pars[2]
  	theta.1.fun   <- pars[3]
  	theta.2.fun   <- pars[4]
  	var.par.1.fun <- pars[5]
  	var.par.2.fun <- pars[6]
  	mu.yr.fun     <- pars[7:length(pars)]	
  }

  # If male or female we keep the maximum fixed so turn off thetas
  # above and turn those on below
  if (theta.const == 1)
  {
  	theta.2.fun   <- (theta.1.fun * (sqrt(1 - cos(2 * pi * max.contr) ^ 2))) /
                     cos(2 * pi * max.contr) 
  }
                   
  # It seems bad but we need to define the bottom bits as global variables
  # so that LogLikelihood can see them
  # Calculate the means given the current update of the parameters
  #print(mean.mnth.coh)
  mean.mnth.coh.fun  <<- matrix(0, nrow = num.months, ncol = no.grps)
   var.mnth.coh.fun  <<- matrix(0, nrow = num.months, ncol = no.grps)
  for (i in seq(1, no.grps))
  {
  	mean.mnth.coh.fun[, i] <<- sapply(months.lst, MeanLength, k0 = k0.fun, theta.1 = theta.1.fun, 
                      	         theta.2 = theta.2.fun, linf = linf.fun , mu.yr = mu.yr.fun,
                      	         yrs.old = yrs.old.par[i], str.mnth = str.mnth.par)       
  }                  
  #print(mean.mnth.coh.fun)
  # Calculate the variances given the current update of the parameters
  for (i in seq(1, no.grps))
  {           
  	var.mnth.coh.fun[, i] <<- sapply(mean.mnth.coh.fun[, i], BscVar, var.par.1 = var.par.1.fun,
                                var.par.2 = var.par.2.fun)
  }
  #print(var.mnth.coh.fun)                
  # Calculate the log likelihood
          
  log.like.full <- sum(sapply(num.months.seq, LogLikelihood1))

  #print(log.like.full)
  # Return the negative of the log likelihood
  
  return(-log.like.full)
}


