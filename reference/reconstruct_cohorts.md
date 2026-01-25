# reconstruct_split_cohort

`reconstruct_split_cohorts()` is used to reconstruct a
\*non-overlapping\* cohort sampled over two years, often with a period
of of zero abundance (e.g. due to egg or adult stages).

`fit_with_offset()` is an internal function used in
\[reconstruct_split_cohort()\] to fit growth functions to find the
optimal cohort offset in a split cohort.

`plot_cohort_fit()` is used to visualize remapped cohorts from
\[reconstruct_split_cohort()\].

## Usage

``` r
reconstruct_split_cohort(
  df = NULL,
  timeCol = "dateID",
  massCol = "massValue",
  massDropThresh = 0.6,
  tStart = 5,
  models = c("vbg", "gompertz", "logistic", "richards"),
  offsetBounds = c(10, 150),
  fallbackGrid = TRUE
)

fit_with_offset(
  dfOrdered = NULL,
  offset = NULL,
  models = c("vbg", "gompertz", "logistic", "richards"),
  tStart = 5
)

plot_cohort_fit(remappedCohort, models = "ensemble", labelPoints = TRUE)
```

## Arguments

- df:

  a data.frame of the sample-level observed masses and densities

- timeCol:

  character. String of the column name depicting the date of each
  sampling event.

- massCol:

  character. String of the column name containing mass data.

- massDropThresh:

  This is the proportional drop in mass used to detect the end and
  beginning of two cohort portions.

- tStart:

  data frame of date information with external predictors for each
  month. There should be a column name identical to all variables in the
  growth equation found in taxaInfo data.frame.

- models:

  character. String vector of the names of models to fit to \\M_t\\. See
  \`details\` for more information.

- offsetBounds:

  integer vector of length = 2. The lower and upper bounds of the offset
  to test for fit

- fallbackGrid:

  logical. If TRUE (default) a grid search procedure will be used if
  \`optim()\` fails.

- dfOrdered:

  the reordered sampling data set from \[reconstruct_split_cohort()\]
  processes

- offset:

  the offset (in days) between the final sample and first sample in two
  cohort portions to be joined

- remappedCohort:

  the reordered object returned from \[reconstruct_split_cohort()\]

- labelPoints:

  logical. Should the points be labelled with the sampling date
  information

## Value

`resonstruct_split_cohort()` returns a dataframe with the original
sampling date, estimated cohort ages, and remapped cohort ages to fit
production methods and remap back to original timescale.

`fit_with_offset()` returns a list of three (3) objects:

fits: the growth model fits to reordered sampling dates

aiccs: the Akaike Information Criterion for each model. This is used to
build ensemble model estimates for remapped cohort 'times'

df_pseudo: a data.frame of the original date, standardized cohort
'times', and t0 corrected ages estimated from the growth models.

`plot_cohort_fit()` returns a ggplot object of the remapped cohort and
growth model fits.

## Details

The initiation of a sampling program can often start when the growth of
a population cohort is already underway, leading to observing only
large, late instars individuals during early sampling events. For
univoltine species, this may be followed by a period of zero abundance
(e.g., due to egg or adult stages) and subsequent sampling of early
instars. To make inferences on cohort parameters and use some
model-based procedures to estimate production, it is useful to
're-order' the sampling dates to create a continuous cohort progression.
This makes some assumptions. For more details on these, see the
\_sampling simulation vignette\_. This function reorders partial cohorts
(e.g., late instars at the beginning and early instars in later sampling
in an annual regime. To do this, the model fits growth functions to
remap split cohorts to a continuous progression in time. The difficult
aspect of this is estimating the relative cohort age of different
sampling events for correct ordering when there may be a long period
between. Crucially, we don't observe the beginning or end of a full
cohort, we only have information from our distinct sampling events. This
function fits growth models to estimate the cohort ages of each sample
and remap into cohort progression.

Models are fit to mean cohort mass at each sampling event, \\M\_{t}\\
and return estimates for asymptotic mass, \\M\_{\infty}\\, time at mass
= 0, \\t_0\\ (or time at growth inflection, \\t\*\\, depending on model
which is used to estimate \\t_0\\), and growth rate, \\k\\. These
Current models available are:

"vbg" The 'regular' von Bertalanffy growth model: \$\$M\_{t} =
M\_{\infty} (1 -e^{-k \cdot (t - t_0)})\$\$ "gompertz" The Gompertz
growth model: \$\$M_t = M\_{\infty} e^{e^{-k \cdot (t - t\*)}}\$\$
"logistic" The logistic growth model: \$\$M_t = M\_{\infty} / (1 + e^{-k
\cdot (t - t\*)})\$\$ "richards" The Richards growth model: \$\$M_t =
M\_{\infty} (1 + 1/D \cdot e^{-k \cdot (t - t\*)})^{-D}\$\$
