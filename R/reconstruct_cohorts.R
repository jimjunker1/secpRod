#' @rdname reconstruct_cohorts
#' @title reconstruct_split_cohort
#' @description
#' \code{reconstruct_split_cohorts()} is used to reconstruct a *non-overlapping* cohort sampled over two years, often with a period of
#' of zero abundance (e.g. due to egg or adult stages).
#' @details
#' The initiation of a sampling program can often start when the growth of a population cohort is already underway, leading to observing only large, late instars individuals during early sampling events.
#' For univoltine species, this may be followed by a period of zero abundance (e.g., due to egg or adult stages) and subsequent sampling of early instars. To make inferences on cohort parameters
#' and use some model-based procedures to estimate production, it is useful to 're-order' the sampling dates to create a continuous cohort progression. This makes some assumptions. For more details on these,
#' see the _sampling simulation vignette_. This function reorders partial cohorts (e.g., late
#' instars at the beginning and early instars in later sampling in an annual regime. To do this, the model fits growth functions to remap split cohorts to a continuous progression in time. The difficult
#' aspect of this is estimating the relative cohort age of different sampling events for correct ordering when there may be a long period between. Crucially, we don't observe the beginning or end of a full
#' cohort, we only have information from our distinct sampling events. This function fits growth models to estimate the cohort ages of each sample and remap into cohort progression.
#'
#' Models are fit to mean cohort mass at each sampling event, \eqn{M_{t}} and return estimates for asymptotic mass, \eqn{M_{\infty}}, time at mass = 0, \eqn{t_0} (or time at growth inflection, \eqn{t*},
#' depending on model which is used to estimate \eqn{t_0}), and growth rate, \eqn{k}. These  Current models available are:
#'
#'      "vbg" The 'regular' von Bertalanffy growth model: \deqn{M_{t} = M_{\infty} (1 -e^{-k \cdot (t - t_0)})}
#'      "gompertz" The Gompertz growth model: \deqn{M_t = M_{\infty} e^{e^{-k \cdot (t - t*)}}}
#'      "logistic" The logistic growth model: \deqn{M_t = M_{\infty} / (1 + e^{-k \cdot (t - t*)})}
#'      "richards" The Richards growth model: \deqn{M_t = M_{\infty} (1 + 1/D \cdot e^{-k \cdot (t - t*)})^{-D}}
#'
#' @param df a data.frame of the sample-level observed masses and densities
#' @param timeCol character. String of the column name depicting the date of each sampling event.
#' @param massCol character. String of the column name containing mass data.
#' @param massDropThresh This is the proportional drop in mass used to detect the end and beginning of two cohort portions.
#' @param tStart This is the guess for the starting value of age (in days) of the smallest mean size sampled.
#' @param models character. String vector of the names of models to fit to estimate \eqn{t_0}. See `details` for more information.
#' @param offsetBounds integer vector of length = 2. The lower and upper bounds of the offset to test for fit
#' @param fallbackGrid logical. If TRUE (default) a grid search procedure will be used if `optim()` fails.
#' @returns \code{resonstruct_split_cohort()} returns a dataframe with the original sampling date, estimated cohort ages, and remapped
#' cohort ages to fit production methods and remap back to original timescale.
#' @import dplyr
#' @importFrom stats coef
#' @importFrom stats optim
#' @export

reconstruct_split_cohort <- function(df = NULL,
                                     timeCol = "dateID",
                                     massCol = "massValue",
                                     massDropThresh = 0.6,
                                     tStart = 5,
                                     models = c("vbg", "gompertz", "logistic", "richards"),
                                     offsetBounds = c(10, 150),
                                     fallbackGrid = TRUE) {
  # declare global vars
  # cohort <- mean_mass <- min_mass <- dateID <- D <- pseudotime <- NULL
  # Step 1a: Aggregate replicate-level masses to daily means
  agg <- df %>%
    group_by(.data[[timeCol]]) %>%
    dplyr::summarise(mean_mass = mean(.data[[massCol]], na.rm = TRUE)) %>%
    ungroup %>%
    dplyr::arrange(.data[[timeCol]])

print('agg')

  t <- agg[[timeCol]]
  W <- agg$mean_mass

  # Step 2a: Detect split cohort via drop in mean mass
  dM <- lead(W) / W
  drop_idx <- which(dM < massDropThresh)[1]
  if (is.na(drop_idx)) drop_idx <- floor(length(t) / 2)  # fallback if no sharp drop

  cohort1 <- agg[1:drop_idx, ]
  cohort2 <- agg[(drop_idx + 1):length(t), ]
print('cohorts')
  # Step 3a: Reorder cohorts so smallest mass cohort comes first
  if (min(cohort1$mean_mass, na.rm = TRUE) > min(cohort2$mean_mass, na.rm = TRUE)) {
    tmp <- cohort1
    cohort1 <- cohort2
    cohort2 <- tmp
  }

  cohort1$cohort <- 1
  cohort2$cohort <- 2
  dfOrdered <- bind_rows(cohort1, cohort2) %>% dplyr::arrange(cohort, .data[[timeCol]])

  cohort_min <-  dfOrdered %>%
    group_by(cohort) %>%
    dplyr::summarise(min_mass = min(mean_mass), .groups = "drop") %>%
    dplyr::arrange(min_mass) %>%
    pull(cohort)

  dfOrdered <-  dfOrdered %>% arrange(cohort, dateID)
print('dfOrdered')
  # 2. Define optimization objective
  obj_fn <- function(offset) {
    out <- fit_with_offset(dfOrdered, offset, models, tStart)
    if (length(out$aiccs) == 0) return(Inf)
    sum(unlist(out$aiccs))
  }
  # debugonce(fit_with_offset)
  # 3. Optimize offset or fallback to grid search
  offset_opt <- tryCatch({
    optim(par = 50, fn = obj_fn, method = "L-BFGS-B",
          lower = offsetBounds[1], upper = offsetBounds[2])$par
  }, error = function(e) {
    # this is currently not operational because the tryCatch will return NULL
    if (fallbackGrid) {
      grid <- seq(offsetBounds[1], offsetBounds[2], by = 1)
      grid_scores <- sapply(grid, obj_fn)
      grid[which.min(grid_scores)]
    } else stop("Optimization failed and grid fallback is off.")
  })
print('optim')
  # 4. Final model fit using optimized offset
  fit_out <- fit_with_offset(dfOrdered, offset_opt, models, tStart)
  fits <- fit_out$fits
  df_pseudo <- fit_out$df_pseudo
print('fit out')
  # 5. Invert model fits to compute pseudotime
  t <- df_pseudo$pseudo_day
  W <- df_pseudo$mean_mass

  cohort_ages <- list()
  weights <- c()

  for (model_name in names(fits)) {
    p <- coef(fits[[model_name]])
    W_cap <- pmin(W, p['Winf'] * 0.999)  # avoid division issues
    t0 <- switch(model_name,
                 vbg = p['t0'],
                 gompertz = p['tStar'] - (1/p['k']) * log(-log(0.0006/p['Winf'])),
                 logistic = p['tStar'] - (1/p['k']) * log(p['Winf']/0.0006 - 1 ),
                 richards = p['tStar'] - (1/p['k']) * log(p['D']*((p['Winf']/0.0006)^(1/p['D']) - 1))
    )
print('t0')
    raw_age <- switch(model_name,
                      vbg = -log(1 - W_cap / p["Winf"]) / p["k"],
                      gompertz = (p["tStar"] - log(-log(W_cap / p["Winf"])) / p["k"]),
                      logistic = (p["tStar"] - log(p["Winf"] / W_cap - 1) / p["k"]),
                      richards = (p["tStar"] - (1 / p['k']) * log(p['D'] * ((p['Winf'] / W_cap)^(1 / p['D']) -1)))
    )
print('raw_age')
    # align such that first pseudotime = abs(t0)
    cohort_ages[[model_name]] <- raw_age - min(raw_age) + abs(t0)
    weights[model_name] <- stats::AIC(fits[[model_name]]) + 2 * 3^2 / (length(t) - 3 - 1)
  }

  # 6. Ensemble weights from AICc
  rel_weights <- exp(-0.5 * (weights - min(weights)))
  rel_weights <- rel_weights / sum(rel_weights)

  pseudo_time <- as.vector(do.call(cbind, cohort_ages) %*% rel_weights)
  pseudo_time <- round(pseudo_time)
print('weights')
  df_remap <- df_pseudo %>%
    mutate(pseudotime = pseudo_time) %>%
    arrange(pseudotime)
print('remap')
  return(list(
    df_remap = df_remap,
    offset = offset_opt,
    weights = rel_weights,
    fits = fits
  ))
}

#' @rdname reconstruct_cohorts
#' @description
#' \code{fit_with_offset()} is an internal function used in [reconstruct_split_cohort()] to fit growth functions to find the optimal cohort offset in a split cohort.
#' @title fit_with_offset
#' @param dfOrdered the reordered sampling data set from [reconstruct_split_cohort()] processes
#' @param offset the offset (in days) between the final sample and first sample in two cohort portions to be joined
#' @param models character. String vector of the names of models to fit to estimate \eqn{t_0}. See `details` for more information.
#' @param tStart data frame of date information with external predictors for each month. There should be a column name identical to all variables in the growth equation found in taxaInfo data.frame.
#' @returns \code{fit_with_offset()} returns a list of three (3) objects:
#' @returns fits: the growth model fits to reordered sampling dates
#' @returns aiccs: the Akaike Information Criterion for each model. This is used to build ensemble model estimates for remapped cohort 'times'
#' @returns df_pseudo: a data.frame of the original date, standardized cohort 'times', and t0 corrected ages estimated from the growth models.
#' @import dplyr
#' @importFrom stats AIC
#' @importFrom stats nls
#' @export

fit_with_offset <- function(dfOrdered = NULL, offset = NULL, models = c("vbg", "gompertz", "logistic", "richards"), tStart = 5) {
  stopifnot("cohort" %in% names(dfOrdered))
# declare global variables
  # cohort <- dateID <- pseudo_day <- NULL
  # Split cohorts
  youngest_dates <- dfOrdered %>%
    dplyr::filter(cohort == min(cohort)) %>%
    arrange(dateID)

  oldest_dates <- dfOrdered %>%
    dplyr::filter(cohort == max(cohort)) %>%
    arrange(dateID)

  # 1. Pseudotime for youngest cohort
  delta_young <- diff(youngest_dates$dateID)
  young_pseudotime <- c(0, cumsum(delta_young))

  # 2. Pseudotime for older cohort
  delta_old <- diff(oldest_dates$dateID)
  if (length(delta_old) == 0) {
    old_pseudotime <- offset
  } else {
    old_pseudotime <- c(offset, offset + cumsum(delta_old))
  }

  # Combine
  df_pseudo <- bind_rows(
    dplyr::mutate(youngest_dates, pseudo_day = young_pseudotime),
    dplyr::mutate(oldest_dates, pseudo_day = old_pseudotime)
  ) %>% arrange(pseudo_day)

  t_all <- df_pseudo$pseudo_day
  W_all <- df_pseudo$mean_mass

  fits <- list()
  aiccs <- c()

  # von Bertalanffy
  if ("vbg" %in% models) {
    m <- try(stats::nls(W_all ~ Winf*(1-exp(-k*(t_all - t0))),
                 start=list(Winf=max(W_all),k=0.01,t0=tStart)), silent=TRUE)
    if (inherits(m,"nls")) {fits$vbg<-m; aiccs["vbg"]<-stats::AIC(m)+2*3^2/(length(t_all)-4)}
  }

  # Gompertz
  if ("gompertz" %in% models) {
    m <- try(stats::nls(W_all ~ Winf*exp(-exp(-k*(t_all - tStar))),
                 start=list(Winf=max(W_all),k=0.01,tStar=tStart)), silent=TRUE)
    if (inherits(m,"nls")) {fits$gompertz<-m; aiccs["gompertz"]<-stats::AIC(m)+2*3^2/(length(t_all)-4)}
  }

  # Logistic
  if ("logistic" %in% models) {
    m <- try(stats::nls(W_all ~ Winf / (1 + exp(-k*(t_all - tStar))),
                 start=list(Winf=max(W_all),k=0.03,tStar=tStart)), silent=TRUE)
    if (inherits(m,"nls")) {fits$logistic<-m; aiccs["logistic"]<-stats::AIC(m)+2*3^2/(length(t_all)-4)}
  }

  # Richards
  if ("richards" %in% models) {
    m <- try(stats::nls(W_all ~ Winf * (1 + 1/D * exp(-k * (t_all - tStar)))^-D,
                 start = list(Winf = max(W_all), k = 0.03, tStar = tStart, D = 50)), silent = TRUE)
    if(inherits(m, 'nls')) {fits$richards <- m; aiccs['richards']<-stats::AIC(m) + 2*3^2/(length(t_all)-4)}
  }

  # exponential
  # if ("exp"%in%models) {
  #   m <- try(stats::nls(W_all ~ k * (t_all - t0),
  #     start=list(k=0.01)), silent=TRUE)
  #   if (inherits(m,"nls")) {fits$schnute<-m; aiccs["exp"]<-stats::AIC(m)+2*2^2/(length(t_all)-3)}
  # }

  return(list(fits=fits, aiccs=aiccs, df_pseudo=df_pseudo))
}

#' @rdname reconstruct_cohorts
#' @description
#' \code{plot_cohort_fit()} is an internal function used in [reconstruct_split_cohort()] to fit growth functions to find the optimal cohort offset in a split cohort.
#' @title plot_cohort_fit
#' @param remappedCohort the reordered object returned from [reconstruct_split_cohort()]
#' @param models character. String vector of the names of models to fit to \eqn{M_t}. See `details` for more information.
#' @param labelPoints logical. Should the points be labelled with the sampling date information
#' @returns \code{plot_cohort_fit()} returns a ggplot object of the remapped cohort and growth model fits.
#' @import dplyr
#' @import ggplot2
#' @importFrom stats coef
#' @export
plot_cohort_fit <- function(remappedCohort, models = "ensemble", labelPoints = TRUE) {

  #declare variables
  # pseudotime <- fitted_mass <- observed_mass <- NULL
  df <- remappedCohort$df_remap
  W <- df$mean_mass
  t <- df$pseudotime


  model_names <- names(remappedCohort$fits)
  fits <- remappedCohort$fits
  weights <- remappedCohort$weights

  pred_df <- tibble(
    pseudotime = t,
    observed_mass = W
  )

  if (models == "ensemble") {
    preds <- matrix(NA, nrow = length(t), ncol = length(model_names))
    colnames(preds) <- model_names

    for (i in seq_along(model_names)) {
      p <- coef(fits[[model_names[i]]])
      W_cap <- pmin(W, p['Winf'] * 0.999)
      preds[, i] <- switch(model_names[i],
                           vbg = p["Winf"] * (1 - exp(-p["k"] * (t - p["t0"]))),
                           gompertz = p["Winf"] * exp(-exp(-p["k"] * (t - p["tStar"]))),
                           logistic = p["Winf"] / (1 + exp(-p["k"]*(t - p["tStar"]))),
                           richards = p["Winf"] * (1 + 1 / p["D"] * exp(-p["k"] * (t - p["tStar"])))^p["D"]
      )
    }
    pred_df$fitted_mass <- as.vector(preds %*% weights)
    model_label <- "Model-averaged fit"

  } else {
    if (!models %in% model_names) stop("Model not found in fits.")
    p <- coef(fits[[models]])
    W_cap <- pmin(W, p['Winf'] * 0.999)
    pred_df$fitted_mass <- switch(models,
                                  vbg = p["Winf"] * (1 - exp(-p["k"] * (t - p["t0"]))),
                                  gompertz = p["Winf"] * exp(-exp(-p["k"] * (t - p["tStar"]))),
                                  logistic = p["Winf"] / (1 + exp(-p["k"]*(t - p["tStar"]))),
                                  richards = p["Winf"] * (1 + 1 / p["D"] * exp(-p["k"] * (t - p["tStar"])))^p["D"]
    )
    model_label <- paste("Fitted:", models)
  }

  ggplot(pred_df, aes(x = pseudotime)) +
    geom_line(aes(y = fitted_mass), color = "blue", linewidth = 1) +
    geom_point(aes(y = observed_mass), color = "black", size = 2) +
    coord_cartesian(xlim = c(0,NA), ylim = c(0,NA))+
    labs(y = "Mean mass (mg)", x = "Pseudotime (days)", title = model_label) +
    theme_minimal(base_size = 13) +
    {if (labelPoints) geom_text(aes(label = df$dateID, y = observed_mass), vjust = -1, size = 3)}
}
