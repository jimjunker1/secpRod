library(TropFishR)
source(here::here("ignore/working/load-packages.R"))
gmean = function(...){
  log_x = log(...)
  meanLog_x = mean(log_x)
  exp(meanLog_x)
}

# load the necessary functions
## functions are loaded in package script.
# source("./R/internals.R")
# source("./R/convert_length_to_mass.R")
# source("./R/wrap_dates.R")
# source("./R/sf_prod.sample.R")
# source("./R/prep_boots.R")
# source("./R/calc_production.R")
# source("./R/calc_prod_sf.R")
# source("./R/plot_cohorts.R")
# source("./R/my_boot.comp.R")


# clean simulated data to workable form

taxaInfo <- data.frame(
  taxonID = c("sppA"),
  massForm = c("afdm_mg~(a*lengthClass^b)*percAsh"),
  a = c(0.0025),
  b = c(2.692),
  percAsh = c(0.958),
  method = c("pb"),
  g.a = c(NA),
  growthForm = c("log(g_d) = 1 - 0.25*log(afdm_mg) - "),
  min.cpi = c(335),
  max.cpi = c(365),
  pb = c("runif(min = 3, max = 8)"),
  min.growth = c(0.001),
  notes = c("This is here for your use. No information will be used, but this column will be maintained in some summaries. See *** for more information.")
)

## load the single species data for working
univoltine <- readRDS("./ignore/working/data/varying_univoltine_samp.rds") %>%
  dplyr::mutate(
    taxonID = "sppA",
    repID = gsub("strain", "", strain)
  ) %>%
  dplyr::filter(stage %ni% c("eggs", "adults")) %>%
  dplyr::select(taxonID, dateID = "date", repID, lengthClass = "length", n_m2)

# create taxa list split on taxonID
taxaList <- univoltine %>% junkR::named_group_split(taxonID)

# convert length to mass
taxaSampleListMass = convert_length_to_mass(taxaSampleList = data.frame(flatten(taxaList)), infoCols = c(2:3), taxaInfo = taxaInfo)

# debugonce(plot_cohorts)
plot_cohorts(taxaSampleListMass, param = 'length', massClass = 'afdm_mg')

# debugonce(calc_production);debugonce(calc_prod_pb);debugonce(pb_prod.sample)
# debugonce(cleanAggDf)
tic();x = calc_production(taxaSampleListMass = taxaSampleListMass, infoCols = c(2:3), taxaInfo = taxaInfo, bootNum = 10, taxaSummary = 'short', wrap = TRUE);toc()

# clean WBT data from Junker and Cross 2014

wbtLenFreq <- readRDS("./ignore/data/WBT_mn.rds") %>%
  rename(taxonID = 'taxon',
         repID = 'SAMPLE',
         dateID = 'DATE',
         lengthClass = 'length_mm') %>%
  dplyr::mutate(lengthClass = ifelse(lengthClass == '<1',0.5,as.numeric(lengthClass))) %>%
  na.omit %>%
  named_group_split(taxonID)

wbtTaxaInfo <- readRDS("./ignore/data/WBT_ab.rds") %>% bind_rows(.id = 'taxonID') %>%
  group_by(taxonID) %>%
  pivot_wider(names_from = 'var_name', values_from = 'value') %>%
  rename(percAsh = 'perc_ash') %>%
  dplyr::mutate(massForm = "afdm_mg~(a*lengthClass^b)*percAsh",
                method = "sf",
                g.a = NA,
                growthForm = "log(g_d) = 1 - 0.25*log(afdm_mg) - ",
                min.cpi = 335,
                max.cpi = 365,
                pb = 5.8,
                min.growth = 0.001,
                notes = "This is here for your use. No information will be used, but this column will be maintained in some summaries. See *** for more information.") %>%
  dplyr::select(taxonID,massForm,a,b,percAsh,method,g.a,growthForm,min.cpi,max.cpi,pb,min.growth,notes ) %>%
  dplyr::mutate(percAsh = case_when(is.na(percAsh) | percAsh == 0 ~ 5.6,
                .default = percAsh))
### working to develop a cohort detection methods
WBTtaxaSampleListMass = convert_length_to_mass(taxaSampleList = wbtLenFreq[[5]], infoCols = c(2:3), taxaInfo = wbtTaxaInfo) %>% dplyr::filter(lengthClass <= 6)
# plot the length frequency data
plot_cohorts(WBTtaxaSampleListMass, param = 'length', massClass = 'afdm_mg')
## estimate Linf from N-Length data

x_Linf = WBTtaxaSampleListMass %>%
  dplyr::select(taxonID, dateID, lengthClass, afdm_mg, n_m2) %>%
  dplyr::filter(n_m2 > 0,
                between(lengthClass,2,7)
                ) %>%
  group_by(lengthClass) %>%
  summarise(count = sum(n_m2))

x_Linf %>%
  ggplot()+
  geom_point(aes(x = lengthClass, y = log(count)))


# debugonce(calc_Linf)
Linf = calc_Linf(x_Linf)
hist(Linf$Linf)
hist(Linf$ZK)


# debugonce(predict_missingN)
predict_missingN(x_Linf, lengthVec = seq(0.1,8,by = 0.1))


# spread the length counts into a single vector
x =  WBTtaxaSampleListMass %>%
  dplyr::select(taxonID, dateID, lengthClass, afdm_mg, n_m2) %>%
  dplyr::filter(n_m2 > 0#,
                # lengthClass <=7
                ) %>%
  group_by(dateID) %>%
  uncount(round(n_m2)) %>%
  dplyr::select(-n_m2) %>%
  ungroup %>%
  mutate(dateID = as.Date(dateID))


restructure_cohorts = function(df,...){
  minL = min(df$lengthClass)
  minDATE = df %>% ungroup %>%
    filter(lengthClass == minL) %>%
    select(dateID) %>%
    slice_min(dateID) %>% unique

  minYEAR = year(minDATE$dateID)
  minYDAY = yday(minDATE$dateID)

  # adjust dates to start at first
  # observation of smallest class

  df_adj = df %>%
    ungroup %>%
    mutate(dateID_adj = yday(dateID)-minYDAY) %>%
    mutate(yday_adj = case_when(dateID_adj < 0 ~ dateID_adj+365,
                                  .default = dateID_adj)) %>%
    arrange(yday_adj)

  return(df_adj)
}
# debugonce(restructure_cohorts)
y = restructure_cohorts(x)

y %>% ggplot()+
  geom_histogram(aes(x = lengthClass, after_stat(..density..)), binwidth = 1)+
  facet_wrap(~yday_adj, dir = 'v' )

min(x$lengthClass)
x %>% filter(lengthClass == min(.$lengthClass)) %>% ungroup %>% select(yday) %>% unique() %>% min()










wbtLFQ = lfqCreate(data = x, Lname = 'lengthClass', Dname = 'dateID')
ELEFAN(wbtLFQ)
Bhattacharya(wbtLFQ)
res= lfqRestructure(wbtLFQ, MA = 3, addl.sqrt = FALSE)
res$peaks_mat
res$rcounts

plot(wbtLFQ_bin2, Fname = 'rcounts')

linf_guess <- max(wbtLFQ$midLengths)/0.95
low_par <- list(Linf = 0.8*linf_guess,
                K = 0.001,
                t_anchor = 0,
                C = 0,
                ts = 0)

up_par <- list(Linf = 1.2*linf_guess,
               K = 1,
               t_anchor = 1,
               C = 1,
               ts = 1)

res_SA <- ELEFAN_SA(wbtLFQ, SA_time = 120, SA_temp = 6e5,
                    MA = 3, seasonalised = TRUE, addl.sqrt = FALSE,
                    init_par = list(Linf = linf_guess,
                                    K = 0.5,
                                    t_anchor = 0.5,
                                    C=0.5,
                                    ts = 0.5),
                    low_par = low_par,
                    up_par = up_par)
res_SA$par
res_SA$Rn_max
plot(wbtLFQ_bin2, Fname = "rcounts",date.axis = "modern", ylim=c(0,10))
lt <- lfqFitCurves(wbtLFQ, par = res_SA$par,
                   draw = TRUE, col = "darkblue", lty = 1, lwd=1.5)
## development of methods to disentangle cohorts
# set initial values
# Gompertz model
cpi = 240 # cpi is the estimated age at max length, first guess
K_init = 6.5 # currently this is an estimate from the length-frequency histogram
M0_init = 0.001 # initial length at t = 0
r_init = 0.05 # initial estimate of growth rates

# first estimate the mean and variance of lengths for each date
# assume there is only a single cohort and normally distributed
# lengths
init_means =
  x %>% group_by(dateID) %>% dplyr::summarise(day = unique(julian(dateID)),
                                              n = n(),
                                              mu = mean(lengthClass,na.rm = TRUE),
                                              sigma = sd(lengthClass, na.rm = TRUE))

# predict growth rates from all sequential sampling dates
# solve for t based on initial guesses for each
gompertzRGR = function(cpi = NULL, K = K_init, M0 = M0_init,...){
  r = -log((log(0.9999999999999999*K/K))/log(M0/K))/cpi
  return(r)
}
gompertzRGR(cpi = 240)

gompertzAge = function(m = NULL, K = K_init, r = gompertzRGR(cpi = 240), M0 = M0_init,...){
  if(any(m == K)){
    m[m==K]<-0.999999999999999*K
  }
  t = -log(log(m/K)/log(M0/K))/r
  return(t)
}
gompertzAge(init_means$mu)

min(gompertzAge(x %>% filter(dateID == as.Date("2008-08-13")) %>% ungroup  %>% select(lengthClass) %>%  unlist, r = gompertzRGR(cpi = 240), M0 = M0_init, K = K_init))
max(gompertzAge(x %>% filter(dateID == as.Date("2008-08-13"), lengthClass <= 6) %>% ungroup  %>% select(lengthClass) %>%  unlist, r = gompertzRGR(cpi = 240), M0 = M0_init, K = K_init))

gompertzAge(m = 6.4999999999,r = gompertzRGR(cpi = 240), M0 = M0_init, K = K_init )

# spread the length counts into a single vector
x =  WBTtaxaSampleListMass %>%
  dplyr::select(taxonID, dateID, lengthClass, afdm_mg, n_m2) %>%
  dplyr::filter(n_m2 > 0,
                lengthClass <=7) %>%
  group_by(dateID) %>%
  uncount(round(n_m2)) %>%
  dplyr::select(-n_m2) %>%
  named_group_split(dateID)

## create a function to estimate the components for each date
# estimate the number of components for each

dateBoots = function(df, lengthCol = 'lengthClass', seed = 1312,...){
  set.seed(seed)
  y = as.vector(unlist(df[lengthCol]))
  lengthClass = makemultdata(y, cuts = quantile(y, (1:9)/10))
  bootybootybooty = my_boot.comp(lengthClass$y, max.comp = 3, B= 1e3, mix.type = "normalmix", sig = 0.095)
  d = bootybootybooty$components
  modes_length = c()
  tryCatch({
  if(d == 1){
    modes_length$mu = mean(y, na.rm = TRUE)
    modes_length$var = sd(y, na.rm = TRUE)
    modes_length$lambda = 1
  } else if(d == 2){
    lambda = 0.5
    mu = quantile(y, c(0.25,0.75))
    sigsqrd = rep(sd(y)/2,d)
    # modes_length = normalmixEM2comp(lengthClass$x, lambda = lambda, mu = mu, sigsqrd = sigsqrd )
    mix = normalmixEM(lengthClass$x, k = 2, arbmean = TRUE, fast = TRUE, maxit = 3e4, maxrestarts = 40)
    modes_length$mu = mix$mu
    modes_length$var = mix$sigma
    modes_length$lambda = mix$lambda
  } else if(d > 2){
    mix = normalmixEM(lengthClass$x, k = d, arbmean = TRUE, arbvar = TRUE)
    modes_length$mu = mix$mu
    modes_length$var = mix$sigma
    modes_length$lambda = mix$lambda
  }}, error = function(e){
    mix = normalmixEM(lengthClass$x, k = d, ECM = TRUE, maxit = 3e4, maxrestarts = 100, epsilon = 1e-03)
    modes_length$mu = mix$mu
    modes_length$var = mix$sigma
    modes_length$lambda = mix$lambda
  }, finally = print("Too many tries"))
  return(list(d = d,
              modes_length = modes_length))
}

# debug(dateBoots)
xDates = x %>% purrr::map(~dateBoots(.x))

## development of methods to disentangle cohorts
# set initial values
# Gompertz model
cpi = 350 # cpi is the estimated age at max length, first guess
linf = 6.5 # currently this is an estimate from the length-frequency histogram
M_0 = 0.01 # initial length at t = 0
r = 0.05 # initial estimate of growth rates

# first estimate the mean and variance of lengths for each date
# assume there is only a single cohort and normally distributed
# lengths
init_means =
  x %>% group_by(dateID) %>% dplyr::summarise(n = n(),
                                              mu = mean(lengthClass,na.rm = TRUE),
                                              sigma = sd(lengthClass, na.rm = TRUE))

# debugonce(calc_production);debugonce(cleanAggDf)
tic();calc_production(taxaSampleListMass = WBTtaxaSampleListMass, infoCols = c(2:3), taxaInfo = wbtTaxaInfo , bootNum = 1e2, taxaSummary = 'short', wrap = TRUE);toc()





source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_log_likelihood_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_mean_length_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_mean_var_optim_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_plot_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_seas_integral_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_seas_root_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_var_ricker_func.R"))
source(here::here("ignore/readings for len-freq function/Zhou et al. supp/mtdgm_bsc/mtdgm_pi_calc_func.R"))
#### mixture models coupled to a von Bertalanffy growth model ####
##### estimate k, Linf, and days to Linf as CPI estimate
# debugonce(convert_length_to_mass)


# Implement Lloyd-Jones framework
## Next steps...need to tweak these functions to allow for daily estimations and convert to monthly aggregates

## Need to update plotting function to automate the date detection for applying this more broadly

## Need to automate the model initialization process to apply more broadly.

## allow the input of other variables than season, i.e., temperature, light, chlorophyll a, etc.

## Can we use the season variable to estimate CPI? or maybe some estimate of CPI based on k and Linf?

####

# Initialize the data for the model
lengths = WBTtaxaSampleListMass %>%
  group_by(taxonID, dateID, lengthClass, afdm_mg) %>%
  dplyr::summarise(n_m2 = sum(n_m2, na.rm = TRUE)) %>%
  dplyr::select(-afdm_mg) %>%
  dplyr::filter(n_m2 > 0) %>%
  uncount(round(n_m2)) %>%
  ungroup %>%
  dplyr::select(lengthClass) %>%
  unlist %>% unname

lfd.year = WBTtaxaSampleListMass %>%
  group_by(taxonID, dateID, lengthClass, afdm_mg) %>%
  dplyr::summarise(n_m2 = sum(n_m2, na.rm = TRUE)) %>%
  dplyr::select(-afdm_mg) %>%
  dplyr::filter(n_m2 > 0) %>%
  uncount(round(n_m2)) %>%
  ungroup %>%
  dplyr::select(dateID) %>%
  dplyr::mutate(year = format(dateID, '%Y')) %>%
  dplyr::select(year) %>% unlist %>% unname

lfd.months = WBTtaxaSampleListMass %>%
  group_by(taxonID, dateID, lengthClass, afdm_mg) %>%
  dplyr::summarise(n_m2 = sum(n_m2, na.rm = TRUE)) %>%
  dplyr::select(-afdm_mg) %>%
  dplyr::filter(n_m2 > 0) %>%
  uncount(round(n_m2)) %>%
  ungroup %>%
  dplyr::select(dateID) %>%
  dplyr::mutate(month = format(dateID, '%m')) %>%
  dplyr::select(month) %>% unlist %>% unname

lfd.days = WBTtaxaSampleListMass %>%
  group_by(taxonID, dateID, lengthClass, afdm_mg) %>%
  dplyr::summarise(n_m2 = sum(n_m2, na.rm = TRUE)) %>%
  dplyr::select(-afdm_mg) %>%
  dplyr::filter(n_m2 > 0) %>%
  uncount(round(n_m2)) %>%
  ungroup %>%
  dplyr::select(dateID) %>%
  dplyr::mutate(yday = yday(dateID)) %>%
  dplyr::select(yday) %>% unlist %>% unname

str.days.yr1 <- min(lfd.days)
end.days.yr1 <- 365
str.days.yr2 <- 0
end.days.yr2 <- 174

lfd.08 <- which((lfd.year == "2008") & (as.numeric(lfd.days) %in% (str.days.yr1:end.days.yr1)))
lfd.09 <- which((lfd.year == "2009") & (as.numeric(lfd.days) %in% (str.days.yr2:end.days.yr2)))

# ------------------------------------------------------------------------------
# Initialise the data for the model
# ------------------------------------------------------------------------------
months.08       <- as.numeric(lfd.months[lfd.08]) - 1 # Jan = 0th month
months.09       <- as.numeric(lfd.months[lfd.09]) + 11
months          <- c(months.08, months.09)
num.months      <- length(names(table(months)))
months.lst      <- as.numeric(names(table(months)))
no.grps.par     <- 2
mu.init.par     <- c(0.05, 6, 1, 7)
# c(K0, Linf, mu1, mu2)
theta.const     <- 0
thetas.par      <- c(-0.2, 0.2)
#thetas.par      <- c(0.68338122,   0.05821895)
var.init.par    <- c(1.5, 0.015)
yrs.old.par     <- c(0, 1)
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
# debugonce(MeanLength)
# debugonce(BscVar)
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
# debugonce(BscPlotNew)
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
#####
x = WBTtaxaSampleListMass %>% dplyr::filter(dateID == as.Date("2009-04-07")) %>%
  dplyr::select(taxonID, lengthClass, afdm_mg, n_m2) %>%
  dplyr::filter(n_m2 > 0,
                lengthClass <=6) %>%
  uncount(round(n_m2)) %>%
  dplyr::select(-n_m2)

mix_length <- mixture(poisson, nmix = 2)
# mix_mass <- mixture(gaussian)

form_length <- bf(lengthClass ~ 1, family = mix_length)
# form_mass <- bf(afdm_mg ~ 1, family = mix_mass)
prior_length = get_prior(form_length, data = x)
# prior_mass = get_prior(form_mass,data = x)

prior_length
# prior_mass

modes_length = normalmixEM(x$lengthClass, k = 2, arbmean = TRUE, arbvar = TRUE)
modes_length$lambda
modes_length$mu
modes_length$sigma
modes_length$loglik


# modes_mass = normalmixEM(x$afdm_mg, k = 2)

lengthPrior_mod = c(prior_string(paste0("normal(",eval(modes_length$mu[1]),", ",eval(modes_length$sigma[1]/3),")"), class = "Intercept", dpar = "mu1"),
              prior_string(paste0("normal(",eval(modes_length$mu[2]),", ",eval(modes_length$sigma[2]/3),")"), class = "Intercept", dpar = "mu2"))

# massPrior_mod = c(prior_string(paste0("normal(",eval(modes_mass$mu[1]),", ",eval(modes_mass$sigma[1]),")"), class = "Intercept", dpar = "mu1"),
              # prior_string(paste0("normal(",eval(modes_mass$mu[2]),", ",eval(modes_mass$sigma[2]),")"), class = "Intercept", dpar = "mu2"))

make_stancode(form_length, family = mix_length, data = x, prior = lengthPrior_mod )

# make_stancode(form_mass, family = mix_mass, data = x, prior = massPrior_mod )

LMmix = brm(form_length,
            data = x,
            family = mix_length,
            prior = lengthPrior_mod,
            chain = 1,
            iter = 1e3)

pp_check(LMmix)
# plot(LMmix)
ppmix = pp_mixture(LMmix, newdata = (x))
hist(ppmix[,1,1])
hist(ppmix[,1,2])
x$group1 = ppmix[,1,1]
x$group2 = ppmix[,1,2]
ggplot(x, aes(x = lengthClass, y = group1))+geom_point()
ggplot(x, aes(x = lengthClass, y = group2))+geom_point()


#### end cohort detection work ####
# debugonce(plot_cohorts)

wbtData = list(sampleInfo = wbtLenFreq,
               taxaInfo = wbtTaxaInfo)



wbtTaxaSampleListMass = purrr::map(wbtLenFreq, ~convert_length_to_mass(.x, infoCols = c(2:3), taxaInfo = wbtTaxaInfo))

purrr::walk(wbtTaxaSampleListMass, ~plot_cohorts(.x, param = 'length', massClass = 'afdm_mg'))

wbtTaxaSampleDfMass = wbttaxaSampleListMass %>% bind_rows()

wbtTaxaInfoList = wbtTaxaInfo %>% named_group_split(taxonID)

tic();wbtTaxaProdList = purrr::map2(wbtTaxaSampleListMass, wbtTaxaInfoList, ~calc_production(taxaSampleListMass = .x, infoCols = c(2:3), taxaInfo = .y, bootNum = 100, taxaSummary = 'short', wrap = TRUE));toc()


predict_missingN = function(df,lengthVec = NULL,...){
  l_nLM = lm(count~lengthClass, data = df)
  pred = predict(l_nLM, newdata = data.frame(lengthClass = lengthVec))
  data.frame(lengthClass = lengthVec,
             count = pred)
}


# # now try a multivoltine species
# taxaInfo <- data.frame(
#   taxonID = c("sppA"),
#   massForm = c("afdm_mg~(a*lengthClass^b)*percAsh"),
#   a = c(0.0025),
#   b = c(2.692),
#   percAsh = c(0.958),
#   method = c("sf"),
#   g.a = c(NA),
#   growthForm = c("log(g_d) = 1 - 0.25*log(afdm_mg) - "),
#   min.cpi = c(335),
#   max.cpi = c(365),
#   pb = c(5.8),
#   min.growth = c(0.001),
#   notes = c("This is here for your use. No information will be used, but this column will be maintained in some summaries. See *** for more information.")
# )
#
# ## load the single species data for working
# multivoltine <- readRDS("./ignore/working/data/multivoltine_samp.rds") %>%
#   dplyr::mutate(taxonID = "sppA") %>%
#   dplyr::select(-doy) %>%
#   group_by(taxonID, date) %>%
#   pivot_longer(cols = c(-taxonID, -date), names_to = 'length', values_to = 'n_m2') %>%
#   dplyr::mutate(length = gsub("multivoltine.","",length)) %>%
#   dplyr::filter(!grepl("eggs|adults", length)) %>%
#   dplyr::mutate(length = as.numeric(gsub("mm","",length))) %>%
#   dplyr::select(taxonID, dateID = "date", lengthClass = "length", n_m2)
#
# # create taxa list split on taxonID
# taxaList <- multivoltine %>% junkR::named_group_split(taxonID)
#
# # convert length to mass
# taxaSampleListMass = convert_length_to_mass(taxaSampleList = data.frame(flatten(taxaList)), infoCols = c(2:3), taxaInfo = taxaInfo)
#
# # debugonce(plot_cohorts)
# plot_cohorts(taxaSampleListMass, param = 'mass', massClass = 'afdm_mg')
