# Performs checks on the estimated divergence times.

This function performs two different checks:

1.  Checks where ancestral nodes without TimeTree information occur most
    often along the tree path: (a) near the root, (b) in the middle,
    or (c) near the tips.

2.  Checks the accuracy of different models for extrapolating the
    TimeTree information, by "erasing" some known points and comparing
    the estimates with the true values.

## Usage

``` r
checkDivergenceTimeEstimates(
  tree,
  timelines = NULL,
  cladeName = NULL,
  models = c("t_gamdec", "t_lin", "t_scam", "t_exp", "t_gam", "t_iso", "t_power")
)
```

## Arguments

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

- timelines:

  An optional named list of data frames, where each name corresponds to
  an NCBI ID. Each element in the list is a data frame containing the
  timeline from the specified taxon to the last common ancestor of all
  cellular organisms. Each row of the dataframe represents one ancestral
  node. This list is computed by this function and returned in
  `timeTreeTimelines`, but because computation can take time, a
  previously computed list can also be given as an optional argument.

- cladeName:

  An optional parameter with the name of the analyzed clade.

- models:

  An optional parameter defining which models will be used for
  extrapolating the data from TimeTree. Models must be listed in order
  of preference. The first model is preferred for setting the divergence
  time. If there is insufficient information to use a model, the next
  model is tried, and so on. Available options: `t_gam` (Generalized
  Additive Model),`t_scam` (Shape Constrained Additive Model),
  `t_gamdec` (Generalized Additive Model with a decreasing monotonicity
  constraint), `t_iso` (Isotonic regression), `t_lin` (Linear
  interpolation), `t_exp` (Exponential model),`t_power` (Power-law
  model).

## Value

A list with two items: 1. `estimatesRes`: the output of the function
`inferDivergenceTimes` (see its documentation for more details); 2.
`accuracySummary`: a dataframe with one row per tested model summarizing
the accuracy of each tested model (e.g., mean, median, maximum, and
minimum error).

## Author

Priscila Biller
