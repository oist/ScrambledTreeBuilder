# Divergence time based on extrapolation

This method infers the divergence time for nodes lacking them based on
the divergence times from other nodes in the tree.

## Usage

``` r
inferDivTimesFromExtrapolation(
  pathToRoot,
  cladeName = NULL,
  verbose = FALSE,
  models = c("t_gamdec", "t_lin", "t_scam", "t_exp", "t_gam", "t_iso", "t_power")
)
```

## Arguments

- pathToRoot:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object containing all nodes along the path from a tip to the root of
  the tree.

- cladeName:

  An optional parameter that describes the dataset.

- verbose:

  Boolean flag (optional; default value is `FALSE`). If `TRUE`, the
  function prints messages about the computation and summary statistics
  for the dataset.

- models:

  An optional parameter defining which models will be used for
  extrapolating the data from TimeTree. Notice that only models with
  sufficient data are included in the results; the amount of information
  required varies by model. Available options: `t_gam` (Generalized
  Additive Model),`t_scam` (Shape Constrained Additive Model),
  `t_gamdec` (Generalized Additive Model with a decreasing monotonicity
  constraint), `t_iso` (Isotonic regression), `t_lin` (Linear
  interpolation), `t_exp` (Exponential model),`t_power` (Power-law
  model).

## Value

A named list with the following: 1. `models`: a named list of the models
actually used to extrapolate the data from TimeTree. A model passed via
the input parameter `models` is omitted if there is insufficient data to
use it. `NULL` is returned if no model was used. 2. `path`: the
[`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
object received as input with divergence-time estimates assigned to
ancestral nodes that lacked them. 3. `known`: nodes used as the basis
for inference. 4. `unknown`: nodes whose divergence times were inferred.
3. `x_grid`: Extra points spanning the region of known and unknown
points where the models were evaluated. 4. `prediction_grids`: A named
list with one element per model, containing predictions for the unknown
points (`unknown`) and the extra points (`x_grid`).

## Details

The method first finds a function describing the relationship between
cumulative branch length and divergence time. It then uses this function
to estimate divergence times for nodes that lack them, given their
cumulative branch lengths.

## Note

See functions `model_gam` (option `t_gam`), `model_scam` (option
`t_scam`), `model_gamdec` (`t_gamdec`), `model_iso` (option `t_iso`),
`model_lin` (option `t_lin`), `model_exp` (option `t_exp`),
`model_power` (option `t_power`) for more details on each model.

As the relationship between cumulative branch length and divergence time
is not necessarily linear and varies for each clade, the model that
seemed to best fit different cases is the Generalized Additive Model
(GAM) with a decreasing monotonicity constraint (`t_gamdec`).

## Author

Priscila Biller
