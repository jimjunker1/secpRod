# Function to visualise the data and the model within loop
# --------------------------------------------------------

BscPlotNew <- function(pars) {
  # Plots the data and the mixture fits as well as the
  # the seasonal function, variance func, linf distribution
  # Args:
  #  pars: Vector containing the parameters to optimise
  #    k0:        Mean  k paramter for VB model
  #    theta.1:   Seasonality parameter 1
  #    theta.2:   Seasonality parameter 1
  #    linf:      Asymptotic length
  #    mu.yr.1:   First month's average length yr 1
  #    mu.yr.2:   First month's average length yr 2T0,T1,LINF,M01,M02
  #    var.pars.x The variance parameters
  #    constrain: flag whether to use a constrained seasonal curve
  # Returns:
  #  A plot of the data and the mixture fits as well as the

  # Declare each of the parameters to names unique to inside the function
  library(ggplot2)
  library(grid)
  library(reshape2)
  constrain  <- pars[1]
  if (constrain == 1)
  {
  	k0         <- pars[2]
  	linf       <- pars[3]
  	theta.1    <- pars[4]
  	var.par.1  <- pars[5]
  	var.par.2  <- pars[6]
  	mu.yr      <- pars[7:length(pars)]
  } else
  {
  	k0         <- pars[2]
  	linf       <- pars[3]
  	theta.1    <- pars[4]
  	theta.2    <- pars[5]
  	var.par.1  <- pars[6]
  	var.par.2  <- pars[7]
  	mu.yr      <- pars[8:length(pars)]
  }
  if ( constrain == 1 ) {
  # If male or female we keep the maximum fixed so turn off thetas
  # above and turn those on below
  theta.2   <- (theta.1 * (sqrt(1 - cos(2 * pi * max.contr)^2))) /
                   cos(2 * pi * max.contr)
  }

  # Declare some plotting parameters
  x.vals <- 0 : 10
  x.lab.name <- c("Jul., 2008", "Aug., 2008", "Sep., 2008", "Oct., 2008", "Nov., 2008",
                  "Dec., 2008", "Jan., 2009", "Feb., 2009", "Mar., 2009", "Apr., 2009",
                  "May, 2009", "Jun., 2009", "Jul., 2009")
  # Calculate the means again for the final likelihood update
  mean.yr <- matrix(0, nrow = num.months, ncol = no.grps)
  var.yr  <- matrix(0, nrow = num.months, ncol = no.grps)
  for (i in seq(1, no.grps))
  {
    mean.yr[, i] <- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1,
                      	         theta.2 = theta.2, linf = linf , mu.yr = mu.yr,
                      	         yrs.old = yrs.old.par[i], str.mnth = str.mnth.par)
  }
  # Calculate the variances again for the final likelihood update
  for (i in seq(1, no.grps))
  {
    var.yr[, i]  <- sapply(mean.yr[, i], BscVar, var.par.1 = var.par.1,
                          var.par.2 = var.par.2)
  }
  pi.all   <- sapply(num.months.seq, PiCalc)
  # Set up the plotting parameters
  dev.new(width = 15, height = 17/1.618, units = "cm")
  if (num.months < 18)
  {
    plot.nrow <- 3
    plot.ncol <- 5
  } else {
  	plot.nrow <- 4
    plot.ncol <- 6
  }
  pushViewport(viewport(layout = grid.layout(plot.nrow, plot.ncol)))
  row.seq <- rep(seq(1, plot.nrow), each  = plot.ncol)
  col.seq <- rep(seq(1, plot.ncol), times = plot.nrow)
  # Plot of data and mixture fits
  for (mm in 1 : num.months) {
  	month.lengths <- lengths[which(months == months.lst[mm])]
    lengths.sub.df <- data.frame(month.lengths)
    m <- ggplot(lengths.sub.df, aes(x = month.lengths)) +
         geom_histogram(aes(y = ..density..),
         fill = rgb(0.5, 0.5, 0.5, 0.5),
         binwidth = 1) + ylim(c(0, 0.66))
    den.yr  <- matrix(0, nrow = length(x.vals), ncol = no.grps)
    for (i in seq(1, no.grps))
    {
    	den.yr[, i] <- pi.all[i, mm] * dnorm(x.vals, mean.yr[mm, i], sqrt(var.yr[mm, i]))
	}
    total.mixture <- rowSums(den.yr)
    x.df     <- data.frame(x.vals, total.mixture, den.yr)
    x.df.mlt <- melt(x.df, id = 1)
    mnth.nm  <- months.lst[mm] + 1
    plot.mm  <- m +
                geom_line(data = x.df.mlt,
                          aes(x = x.vals,
                              y = value,
                              group = variable,
                              lty   = variable),
                              size = 0.5) +
                xlab(x.lab.name[mnth.nm]) + theme(legend.position = "none")
    #print(c(row.seq[mm], col.seq[mm]))
    print(plot.mm, vp = viewport(layout.pos.row = row.seq[mm],
          layout.pos.col = col.seq[mm]))
  }
  if (num.months < 18)
  {
    lpr.seas <- 3
    lpr.var  <- 3
    lpr.linf <- 3
    seas.plt.col <- 2
    var.plt.col  <- 3
    var.ds.plt.col  <- 4
  } else {
  	lpr.seas <- 4
    lpr.var  <- 4
    lpr.linf <- 4
    seas.plt.col <- 3
    var.plt.col  <- 4
    var.ds.plt.col  <- 5
  }
  # Plot the seasonal function
  # --------------------------
	x.vals.2  <- seq(0, 1, length.out = 100)
	yy <- k0 + theta.1 * cos(2 * pi * x.vals.2) +
	               theta.2 * sin(2 * pi * x.vals.2)
	seas.func <- ((yy + abs(yy))/2)
	df.seas <- data.frame(x.vals.2, seas.func)
	max.seas <- max(seas.func) + 0.05
	seas.plot <- ggplot(df.seas, aes(x = x.vals.2 , y = seas.func)) +
	               xlab("Fraction of a year since Jan. 1") +
	               ylab(expression(paste("g(", bold(theta), ", t)"))) +
	               ylim(c(0, max.seas))
    print(seas.plot + geom_line(size = 1),
	      vp = viewport(layout.pos.row = lpr.seas,
	                    layout.pos.col = seas.plt.col))
	# Plot the variance funtion
	# -------------------------
	#x.vals <- 0 : 400
	#var.par.1 <- 19
	#var.par.2 <- 0.015
	var.linf.x <-  sapply(x.vals, BscVar, var.par.1 =  var.par.1,
	                     var.par.2 = var.par.2)
	df.var   <- data.frame(x.vals, var.linf.x)
	var.plot <- ggplot(df.var, aes(x = x.vals, y = var.linf.x)) +
	               xlab("Mean length (mm)") + ylab("Variance")
	#var.plot + geom_line(size = 1)
	print(var.plot + geom_line(size = 1), vp = viewport(layout.pos.row = lpr.var,
	        layout.pos.col = var.plt.col))
	# # Plot the distribution of linf
	# # -----------------------------
	# var.linf <- BscVar(var.par.1, var.par.2, linf)
	# var.linf.vals <- rnorm(1000000, linf, sqrt(var.linf))
	# var.linf.dist <- data.frame(var.linf.vals)
	# m <- ggplot(var.linf.dist, aes(x = var.linf.vals)) +
	     # xlab("Length (mm)") +  stat_density(geom = "line", size = 1)
	# print(m,
	      # vp = viewport(layout.pos.row = lpr.linf,
	      # layout.pos.col = var.ds.plt.col))

}



























