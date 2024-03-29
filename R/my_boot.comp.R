#' @title my_boot.comp
#' @param y The raw data for multmix, mvnormalmix, normalmix, and repnormmix and the response values for logisregmix, poisregmix, and regmix. See the documentation concerning their respective EM algorithms for specific structure of the raw data.
#' @param x The predictor values required only for the regression mixtures logisregmix, poisregmix, and regmix. A column of 1s for the intercept term must not be included! See the documentation concerning their respective EM algorithms for specific structure of the predictor values.
#' @param N An n-vector of number of trials for the logistic regression type logisregmix. If NULL, then N is an n-vector of 1s for binary logistic regression.
#' @param max.comp The maximum number of components to test for. The default is 2. This function will perform a test of k-components versus (k+1)-components sequentially until we fail to reject the null hypothesis. This decision rule is governed by the calculated p-value and sig.
#' @param B The number of bootstrap realizations of the likelihood ratio statistic to produce. The default is 100, but ideally, values of 1000 or more would be more acceptable.
#' @param sig The significance level for which to compare the p-value against when performing the test of k-components versus (k+1)-components.
#' @param arbmean If FALSE, then a scale mixture analysis can be performed for mvnormalmix, normalmix, regmix, or repnormmix. The default is TRUE.
#' @param arbvar If FALSE, then a location mixture analysis can be performed for mvnormalmix, normalmix, regmix, or repnormmix. The default is TRUE.
#' @param mix.type The type of mixture analysis you wish to perform. The data inputted for y and x depend on which type of mixture is selected. logisregmix corresponds to a mixture of logistic regressions. multmix corresponds to a mixture of multinomials with data determined by the cut-point method. mvnormalmix corresponds to a mixture of multivariate normals. normalmix corresponds to a mixture of univariate normals. poisregmix corresponds to a mixture of Poisson regressions. regmix corresponds to a mixture of regressions with normal components. regmix.mixed corresponds to a mixture of regressions with random or mixed effects. repnormmix corresponds to a mixture of normals with repeated measurements.
#' @param hist An argument to provide a matrix plot of histograms for the boostrapped likelihood ratio statistic.
#' @param ... Additional arguments passed to the various EM algorithms for the mixture of interest
#' @returns \code{boot.comp} returns a list with items:\item{p.values}{The p-values for each test of k-components versus (k+1)-components.} \item{log.lik}{The B bootstrap realizations of the likelihood ratio statistic.} \item{obs.log.lik}{The observed likelihood ratio statistic for each test which is used in determining the p-values.}
#' @import stats
#' @importFrom mixtools regmixEM
#' @importFrom mixtools repnormmixEM
#' @importFrom mixtools dmvnorm
#' @importFrom mixtools mvnormalmixEM
#' @importFrom mixtools rmvnorm
#' @importFrom mixtools rmvnormmix
#' @importFrom mixtools normalmixEM
#' @importFrom mixtools ldmult
#' @importFrom mixtools multmixEM
#' @importFrom mixtools logisregmixEM
#' @importFrom mixtools poisregmixEM
#' @export
#'
my_boot.comp = function (y, x = NULL, N = NULL, max.comp = 2, B = 100, sig = 0.05,
          arbmean = TRUE, arbvar = TRUE, mix.type = c("logisregmix",
                                                      "multmix", "mvnormalmix", "normalmix", "poisregmix",
                                                      "regmix", "regmix.mixed", "repnormmix"), hist = TRUE,
          ...)
{
  mix.type <- match.arg(mix.type)
  k = max.comp
  p = 0
  sigtest = 1
  Q.star = list()
  i = 0
  if (mix.type == "regmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = lm(y ~ x)
          beta = coef(H0.fit)
          Q0[i] = as.numeric(logLik(H0.fit))
          H1.fit = try(regmixEM(y = y, x = x, k = (i +
                                                     1), arbmean = arbmean, arbvar = arbvar,
                                ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0)
              w = 1
            else w = 2
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            xy.sd = sqrt(sum(H0.fit$res^2)/(length(y) -
                                              2))
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rnorm(length(y), mean = xbeta, sd = xy.sd)
          xy.simout = lm(y.sim ~ x)
          em.out = try(regmixEM(y = y.sim, x = x, k = (i +
                                                         1), arbmean = arbmean, arbvar = arbvar,
                                ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - as.numeric(logLik(xy.simout)))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(regmixEM(y = y, x = x, k = i,
                                arbmean = arbmean, arbvar = arbvar, ...),
                       silent = TRUE)
          H1.fit = try(regmixEM(y = y, x = x, k = (i +
                                                     1), arbmean = arbmean, arbvar = arbvar,
                                ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) {
              scale = H0.fit$scale
              beta = matrix(rep(H0.fit$beta, i), ncol = i)
            }
            else {
              scale = 1
            }
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0)
              w = 1
            else w = 2
          }
          beta.new = H0.fit$beta
          xbeta.new = cbind(1, x) %*% beta.new
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                          mean = xbeta.new, sd = ((scale * H0.fit$sigma)[wt[,
                                                                                                            i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                            mean = xbeta.new[i, (wt[, i] == 1)],
                                                            sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                            mean = xbeta.new[i, (wt[, i] == 1)],
                                                            sd = H0.fit$sigma[wt[, i] == 1]))
            }
          }
          em.out.0 = try(regmixEM(y = y.sim, x = x,
                                  k = i, arbmean = arbmean, arbvar = arbvar,
                                  ...), silent = TRUE)
          em.out.1 = try(regmixEM(y = y.sim, x = x,
                                  k = (i + 1), arbmean = arbmean, arbvar = arbvar,
                                  ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "repnormmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          dens = dnorm(y, mean = mean(y), sd = sd(y))
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(repnormmixEM(x = y, k = (i +
                                                  1), arbmean = arbmean, arbvar = arbvar,
                                    ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rmvnormmix(nrow(y), mu = rep(mean(y),
                                               ncol(y)), sigma = rep(sd(y), ncol(y)))
          dens.sim = dnorm(y.sim, mean = mean(y), sd = sd(y))
          x.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(repnormmixEM(x = y.sim, k = (i +
                                                      1), arbmean = arbmean, arbvar = arbvar,
                                    ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - x.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(repnormmixEM(x = y, k = i, arbmean = arbmean,
                                    arbvar = arbvar, ...), silent = TRUE)
          H1.fit = try(repnormmixEM(x = y, k = (i +
                                                  1), arbmean = arbmean, arbvar = arbvar,
                                    ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE)
              scale = H0.fit$scale
            else scale = 1
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y),
                                                        mean = H0.fit$mu, sd = ((scale * H0.fit$sigma)[wt[,
                                                                                                          i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y),
                                                          mean = H0.fit$mu[wt[, i] == 1], sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:ncol(y), function(i) rnorm(nrow(y),
                                                          mean = H0.fit$mu[wt[, i] == 1], sd = H0.fit$sigma[wt[,
                                                                                                               i] == 1]))
            }
          }
          em.out.0 = try(repnormmixEM(x = y.sim, k = i,
                                      arbmean = arbmean, arbvar = arbvar, ...),
                         silent = TRUE)
          em.out.1 = try(repnormmixEM(x = y.sim, k = (i +
                                                        1), arbmean = arbmean, arbvar = arbvar,
                                      ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "regmix.mixed") {
    if (is.list(y)) {
      if (length(y) != length(x))
        stop("Number of elements in lists for x and y must match!")
    }
    tt = sapply(1:length(x), function(i) x[[i]][, 1])
    beta = t(sapply(1:length(y), function(i) lm(y[[i]] ~
                                                  x[[i]])$coef))
    y = beta
    mix.type = "mvnormalmix"
  }
  if (mix.type == "mvnormalmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          y.mean = apply(y, 2, mean)
          y.cov = cov(y)
          dens = dmvnorm(y, mu = y.mean, sigma = y.cov)
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(mvnormalmixEM(x = y, k = (i +
                                                   1), arbmean = arbmean, arbvar = arbvar,
                                     ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0)
              w = 1
            else w = 2
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rmvnorm(nrow(y), mu = apply(y, 2,
                                              mean), sigma = y.cov)
          dens.sim = dmvnorm(y.sim, mu = y.mean, sigma = y.cov)
          y.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(mvnormalmixEM(x = y.sim, k = (i +
                                                       1), arbmean = arbmean, arbvar = arbvar,
                                     ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - y.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(mvnormalmixEM(x = y, k = i, ...),
                       silent = TRUE)
          H1.fit = try(mvnormalmixEM(x = y, k = (i +
                                                   1), arbmean = arbmean, arbvar = arbvar,
                                     ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE) {
              H0.fit$mu = lapply(1:i, function(l) H0.fit$mu)
            }
            if (arbvar == FALSE) {
              H0.fit$sigma = lapply(1:i, function(l) H0.fit$sigma)
            }
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j <- 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(nrow(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1,
                                                            mu = H0.fit$mu, sigma = H0.fit$sigma[wt[,
                                                                                                    i] == 1][[1]])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1,
                                                              mu = H0.fit$mu[wt[, i] == 1][[1]], sigma = H0.fit$sigma)))
            }
            else {
              y.sim = t(sapply(1:nrow(y), function(i) rmvnorm(1,
                                                              mu = H0.fit$mu[wt[, i] == 1][[1]], sigma = H0.fit$sigma[wt[,
                                                                                                                         i] == 1][[1]])))
            }
          }
          em.out.0 = try(mvnormalmixEM(x = y.sim, k = i,
                                       arbmean = arbmean, arbvar = arbvar, ...),
                         silent = TRUE)
          em.out.1 = try(mvnormalmixEM(x = y.sim, k = (i +
                                                         1), arbmean = arbmean, arbvar = arbvar,
                                       ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "normalmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          dens = dnorm(y, mean = mean(y), sd = sd(y))
          Q0[i] = sum(log(dens[dens > 0]))
          H1.fit = try(normalmixEM(x = y, k = (i + 1),
                                   arbmean = arbmean, arbvar = arbvar, ...),
                       silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rnorm(length(y), mean = mean(y), sd = sd(y))
          dens.sim = dnorm(y.sim, mean = mean(y), sd = sd(y))
          x.simout = sum(log(dens.sim[dens.sim > 0]))
          em.out = try(normalmixEM(x = y.sim, k = (i +
                                                     1), arbmean = arbmean, arbvar = arbvar,
                                   ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - x.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(normalmixEM(x = y, k = i, arbmean = arbmean,
                                   arbvar = arbvar, ...), silent = TRUE)
          H1.fit = try(normalmixEM(x = y, k = (i + 1),
                                   arbmean = arbmean, arbvar = arbvar, ...),
                       silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            if (arbmean == FALSE)
              scale = H0.fit$scale
            else scale = 1
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          if (arbmean == FALSE) {
            y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                          mean = H0.fit$mu, sd = ((scale * H0.fit$sigma)[wt[,
                                                                                                            i] == 1])))
          }
          else {
            if (arbvar == FALSE) {
              y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                            mean = H0.fit$mu[(wt[, i] == 1)], sd = H0.fit$sigma))
            }
            else {
              y.sim = sapply(1:length(y), function(i) rnorm(1,
                                                            mean = H0.fit$mu[(wt[, i] == 1)], sd = H0.fit$sigma[wt[,
                                                                                                                   i] == 1]))
            }
          }
          em.out.0 = try(normalmixEM(x = y.sim, k = i,
                                     arbmean = arbmean, arbvar = arbvar, ...),
                         silent = TRUE)
          em.out.1 = try(normalmixEM(x = y.sim, k = (i +
                                                       1), arbmean = arbmean, arbvar = arbvar,
                                     ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "multmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          m = apply(y, 1, sum)
          n.i = apply(y, 2, sum)
          theta = n.i/sum(n.i)
          Q0[i] = sum(log(exp(apply(y, 1, ldmult, theta = theta))))
          H1.fit = try(multmixEM(y = y, k = (i + 1),
                                 ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = matrix(0, ncol = ncol(y), nrow = nrow(y))
          for (l in 1:length(m)) {
            y.sim[l, ] <- rmultinom(1, size = m[l],
                                    prob = theta)
          }
          theta.sim = apply(y.sim, 2, sum)/sum(apply(y.sim,
                                                     2, sum))
          y.simout = sum(log(exp(apply(y.sim, 1, ldmult,
                                       theta = theta))))
          em.out = try(multmixEM(y = y.sim, k = (i +
                                                   1), ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - y.simout)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(multmixEM(y = y, k = i, ...),
                       silent = TRUE)
          H1.fit = try(multmixEM(y = y, k = (i + 1),
                                 ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            theta = H0.fit$theta
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(nrow(y), size = 1, prob = H0.fit$lambda)
          new.y.sim = t(sapply(1:nrow(y), function(i) rmultinom(1,
                                                                size = n.i, prob = H0.fit$theta[(wt[, i] ==
                                                                                                   1), ])))
          em.out.0 = try(multmixEM(y = new.y.sim, k = i,
                                   ...), silent = TRUE)
          em.out.1 = try(multmixEM(y = new.y.sim, k = (i +
                                                         1), ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            Q.star[[i]][j]
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "logisregmix") {
    if (is.null(N))
      stop("Number of trials must be specified!")
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = glm(cbind(y, N - y) ~ x, family = binomial())
          Q0[i] = logLik(H0.fit)
          H1.fit = try(logisregmixEM(y = y, x = x, N = N,
                                     k = (i + 1), ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rbinom(length(y), size = N, prob = inv.logit(xbeta))
          xy.simout = glm(cbind(y.sim, N - y.sim) ~
                            x, family = binomial())
          em.out = try(logisregmixEM(y = y.sim, x = x,
                                     N = N, k = (i + 1), ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - logLik(xy.simout))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(logisregmixEM(y = y, x = x, N = N,
                                     k = i, ...), silent = TRUE)
          H1.fit = try(logisregmixEM(y = y, x = x, N = N,
                                     k = (i + 1), ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = H0.fit$beta
            xbeta = cbind(1, x) %*% beta
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          y.sim = sapply(1:length(y), function(i) rbinom(1,
                                                         size = N[i], prob = inv.logit(xbeta)[, (wt[,
                                                                                                    i] == 1)]))
          em.out.0 = try(logisregmixEM(y = y.sim, x = x,
                                       N = N, k = i, ...), silent = TRUE)
          em.out.1 = try(logisregmixEM(y = y.sim, x = x,
                                       N = N, k = (i + 1), ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (mix.type == "poisregmix") {
    Q0 = 0
    Q1 = 0
    obs.Q = 0
    i = 1
    while (sigtest == 1 && i <= k) {
      Q.star[[i]] = 0
      if (i == 1) {
        w = 1
        while (w == 1) {
          H0.fit = glm(y ~ x, family = poisson())
          Q0[i] = logLik(H0.fit)
          H1.fit = try(poisregmixEM(y = y, x = x, k = (i +
                                                         1), ...), silent = TRUE)
          if (inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = coef(H0.fit)
            xbeta = cbind(1, x) %*% beta
            j = 0
          }
        }
        while (j < B) {
          j = j + 1
          y.sim = rpois(length(y), lambda = exp(xbeta))
          xy.simout = glm(y.sim ~ x, family = poisson())
          em.out = try(poisregmixEM(y = y.sim, x = x,
                                    k = (i + 1), ...), silent = TRUE)
          if (inherits(em.out, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out$loglik - logLik(xy.simout))
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      else {
        w = 1
        while (w == 1) {
          H0.fit = try(poisregmixEM(y = y, x = x, k = i,
                                    ...), silent = TRUE)
          H1.fit = try(poisregmixEM(y = y, x = x, k = (i +
                                                         1), ...), silent = TRUE)
          if (inherits(H0.fit, "try-error", which = TRUE) ||
              inherits(H1.fit, "try-error", which = TRUE)) {
            w = 1
          }
          else {
            Q0[i] = H0.fit$loglik
            Q1[i] = H1.fit$loglik
            obs.Q[i] = 2 * (Q1[i] - Q0[i])
            if (obs.Q[i] < 0) {
              w = 1
            }
            else {
              w = 2
            }
            beta = H0.fit$beta
            xbeta = cbind(1, x) %*% beta
          }
          j = 0
        }
        while (j < B) {
          j = j + 1
          wt = rmultinom(length(y), size = 1, prob = H0.fit$lambda)
          y.sim = sapply(1:length(y), function(i) rpois(1,
                                                        lambda = exp(xbeta)[, (wt[, i] == 1)]))
          em.out.0 = try(poisregmixEM(y = y.sim, x = x,
                                      k = i, ...), silent = TRUE)
          em.out.1 = try(poisregmixEM(y = y.sim, x = x,
                                      k = (i + 1), ...), silent = TRUE)
          if (inherits(em.out.0, "try-error", which = TRUE) ||
              inherits(em.out.1, "try-error", which = TRUE)) {
            j = j - 1
          }
          else {
            Q.star[[i]][j] = 2 * (em.out.1$loglik -
                                    em.out.0$loglik)
            if (Q.star[[i]][j] < 0) {
              j = j - 1
            }
          }
        }
      }
      p[i] = mean(Q.star[[i]] >= obs.Q[i])
      sigtest = (p[i] < sig)
      i = i + 1
    }
  }
  if (hist) {
    if (length(p) == 2) {
      par(mfrow = c(1, 2))
      for (i in 1:length(p)) {
        hist(Q.star[[i]], xlab = c("Bootstrap Likelihood",
                                   "Ratio Statistic"), main = paste(i, "versus",
                                                                    i + 1, "Components"))
        segments(obs.Q[i], 0, obs.Q[i], B, col = 2,
                 lwd = 2)
      }
    }
    else {
      g = ceiling(sqrt(length(p)))
      par(mfrow = c(g, g))
      for (i in 1:length(p)) {
        hist(Q.star[[i]], xlab = c("Bootstrap Likelihood",
                                   "Ratio Statistic"), main = paste(i, "versus",
                                                                    i + 1, "Components"))
        segments(obs.Q[i], 0, obs.Q[i], B, col = 2,
                 lwd = 2)
      }
    }
  }
  if (p[length(p)] < sig) {
    cat("Decision: Select", length(p) + 1, "Component(s)",
        "\n")
    comps = (length(p) + 1)
  }
  else {
    cat("Decision: Select", length(p), "Component(s)", "\n")
    comps = (length(p))
  }
  list(p.values = p, log.lik = Q.star, obs.log.lik = obs.Q, components = comps)
}
