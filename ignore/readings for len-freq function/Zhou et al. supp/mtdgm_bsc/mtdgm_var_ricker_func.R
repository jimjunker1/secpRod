# Variance function for calcuating the variance for the current month
# -------------------------------------------------------------------


BscVar <- function(var.par.1, var.par.2, mean.length) {
	
  # Computes the variance of the distribution for the current month
  # based on a bespoke varinace function that is positive in x and y.
  #
  # Args:
  #  var.par.1:   Variance parameter 1
  #  var.par.2:   Variance parameter 2
  #  mean.length: Mean length for current month
  #
  # Returns:
  #  Variance for each month
  
  variance = max(var.par.1 * mean.length * exp(-var.par.2 * mean.length), 1)


  return(variance)
}