# Function to calculate roots of the integral of the VB seasonal curve
# ---------------------------------------------------------------------

SeasIntegFunc <- function(k0, theta.1, theta.2, str.time, end.time) {
  # Calculates the canonical seasonal function for the VB model
  #
  # Args:
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  str.time  start time for integral
  #  end.time  end time for integral
  #
  # Returns:
  #  The resultant value of the seasonal function
  
  integral <- k0 * (end.time - str.time) + 
              (theta.1 / (2 * pi)) * (sin(2 * pi * end.time) - sin(2 * pi * str.time)) - 
              (theta.2 / (2 * pi)) * (cos(2 * pi * end.time) - cos(2 * pi * str.time))
              
  return(integral)
}
