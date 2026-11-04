# Estimate divergence time for ancestral nodes

This method infers the divergence time of ancestral nodes based on the
information from TimeTree when available, and by extrapolating this
information to other nodes without this information when possible.

## Usage

``` r
inferDivergenceTimes(
  tree,
  timelines = NULL,
  cladeName = NULL,
  savePlots = FALSE,
  verbose = FALSE,
  models = c("t_gamdec", "t_lin", "t_scam", "t_exp", "t_gam", "t_iso")
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

- savePlots:

  Boolean flag (optional; default value is `FALSE`). If `TRUE`, the
  function saves `.png` and `.svg` files containing several plots (e.g.
  relationship between cumulative branch length and divergence time,
  stats about the nodes for each path from tip to root, etc.)

- verbose:

  Boolean flag (optional; default value is `FALSE`). If `TRUE`, the
  function prints messages about the computation and summary statistics
  for the dataset.

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

A named list with the following information: 1. `tree`: A
[`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
object with the following additional columns: `divergenceTime` (the
estimated divergence time for the node, either coming from TimeTree or
extrapolated from TimeTree information), `isTimeExtrapolated` (`FALSE`
if the divergence time comes directly from TimeTree, `TRUE` if it is
extrapolated from TimeTree), `cumBranchLen` (the sum of the edge lengths
from the root to the node), `isSingleSpecies` (`TRUE` if all tips share
the same `NCBI_id`; false otherwise), plus additional columns, one per
model used, named after the model and containing the divergence times
inferred by that model. 2. `timeTreeTimelines`: A list of data frames,
where each name corresponds to an NCBI ID (from column `NCBI_id`). Each
element in the list is a data frame containing the timeline from the
specified taxon to the last common ancestor of all cellular organisms.
Each row of the dataframe represents one ancestral node. 3.
`timeTreeDivergenceTimes`: A named list of data frames, where each
element corresponds to an ancestral node in the tree for which a
divergence time is available from TimeTree (names correspond to column
`node`). Each element has a dataframe (computed with the function
`inferDivTimesFromTimeTree`), with one row per informative pair of
leaves and the divergence time inferred from that pair. Notice that an
ancestral node may have multiple pairs of leaves that yield different
divergence times, indicating a possible disagreement between the
TimeTree phylogeny and the phylogeny from ScrambledTreeBuilder. 4.
`pathsLeafToRoot`: All paths from tips to root used for the divergence
time estimation, sorted by the amount of information from TimeTree they
have.

## Note

See functions `model_gam` (option `t_gam`), `model_scam` (option
`t_scam`), `model_gamdec` (`t_gamdec`), `model_iso` (option `t_iso`),
`model_lin` (option `t_lin`), `model_exp` (option `t_exp`),
`model_power` (option `t_power`) for more details on each model.

## Author

Priscila Biller
