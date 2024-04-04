# Function to calculate roots of the seasonal function for the VB model
# ---------------------------------------------------------------------

SeasFunc <- function(k0, theta.1, theta.2, time.var) {
  # Calculates the canonical seasonal function for the VB model
  #
  # Args:
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  time.var  time variable
  #
  # Returns:
  #  The resultant value of the seasonal function
  
  return(k0 + theta.1 * cos(2 * pi * time.var) + theta.2 * sin(2 * pi * time.var))
}


SeasRootCalc <- function(k0, theta.1, theta.2, yrs.old) {
	
  # Computes the roots of the seasonal function used 
  # in VB growth models. This allows for easy integration
  # of the 
  #
  # Args:
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  yrs.old:  How old the individuals in this group are
  #
  # Returns:
  #  The two roots between 0 and 1 
  
  # Define a set of aux variables
  
  a = theta.1^2 + theta.2^2
  b = 2 * k0 * theta.1
  c = -theta.2^2 + k0^2
  
  # Find the set of aux variables solutions
  
  u1 = (-b + sqrt(b^2 -4 * a * c)) / (2 * a)
  u2 = (-b - sqrt(b^2 -4 * a * c)) / (2 * a)
  
  # Map them back via arc cos
  
  r11 = 1 - acos(u1) / (2 * pi) + yrs.old
  r12 =     acos(u1) / (2 * pi) + yrs.old
  r21 = 1 - acos(u2) / (2 * pi) + yrs.old
  r22 =     acos(u2) / (2 * pi) + yrs.old
  
  # Find those that satisfy the roots of our function
  
  roots  = c(r11, r12, r21, r22)
  g.root = SeasFunc(k0, theta.1, theta.2, roots)
  g.min  = round(g.root)
  r1     = min(roots[which(g.min == 0)])
  r2     = max(roots[which(g.min == 0)])
  
  # Return the roots
  
  return(c(r1, r2))
}
  