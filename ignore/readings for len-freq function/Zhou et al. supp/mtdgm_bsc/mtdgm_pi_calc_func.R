# Function to calculate pis for each month
# ----------------------------------------
	
PiCalc <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  mm: Index for which month we are up to. Not the month itself
  #
  # Returns:
  #  A set of pis for each group for the month
  
  # Subset the lengths
  #mm = 1
  lengths.sub    <- lengths[which(months == months.lst[mm])]
  
  # Calculate the inclusion probablities 
  top <- matrix(0, nrow = length(lengths.sub), ncol = no.grps)
  bot <- array(0,  length(lengths.sub))
  for (i in seq(1, no.grps))
  {
    top[, i] <- pis[mm, i] * dnorm(lengths.sub, mean.mnth.coh.mn[mm, i], sqrt(var.mnth.coh.mn[mm, i]))	
  }		
  bot  <- rowSums(top)
  # Calculate the tau scores for each group
  taus <- matrix(0, nrow = length(lengths.sub), ncol = no.grps)
  for (i in seq(1, no.grps))
  {
    taus[, i] <- top[, i] / bot	
  }		
  # Return the pis for each group	  	
  pis.mm <- colSums(taus) / dim(taus)[1]
  # Return set for that months
  return(pis.mm)
}