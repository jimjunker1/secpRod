# Function to calculate log likelihood
# ------------------------------------
	
LogLikelihood1 <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  mm: Index for which month we are up to. Not the month itself
  #
  # Returns:
  #  The log likelihood for each month
  #mm = 1
  #print(mean.mnth.coh.fun)
  lengths.sub    <- lengths[which(months == months.lst[mm])]
  like <- 0.0
  #print(mean.mnth.coh)
  for (i in seq(1, no.grps))
  {
  	#i = 1
    like.grp <- pis[mm, i] * dnorm(lengths.sub, mean.mnth.coh.fun[mm, i], 
                sqrt(var.mnth.coh.fun[mm, i]))
    like     <- like + like.grp
  }
  log.like <- sum(log(like))
  return(log.like)

}

LogLikelihood2 <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  mm: Index for which month we are up to. Not the month itself
  #
  # Returns:
  #  The log likelihood for each month
  #mm = 1
  lengths.sub    <- lengths[which(months == months.lst[mm])]
  like <- 0.0
  #print(mean.mnth.coh)
  for (i in seq(1, no.grps))
  {
  	#i = 1
    like.grp <- pis[mm, i] * dnorm(lengths.sub, mean.mnth.coh.mn[mm, i], 
                sqrt(var.mnth.coh.mn[mm, i]))
    like     <- like + like.grp
  }
  log.like <- sum(log(like))
  return(log.like)

}
