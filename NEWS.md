# ScrambledTreeBuilder dev

* Standardise function argument names: use `pairwise_data` for data frames,
  `pairwise_matrix` for matrices, and `species_tree` for trees.  The S7 class
  is named `stbTree` (avoiding conflict with `ape::speciesTree`).
* `MRCA_2D_plot()` now runs `MRCAs()` internally and adds hover information to
  make the plot ready for Plotly (interactive features contributed by Takahiro
  Fujita).
* New `computeENR()` and `cladeENRtable()` functions for detecting outlier
  clades.
* `visualizeTree()` can color the branches of a tree by a species trait
  stored on its tips, such as the chromosome number (`trait` option), and
  print the trait values next to the species names.  New
  `ancestralStates()` function that reconstructs ancestral values of a
  species trait by parsimony or by maximum likelihood (with `ape::ace()`,
  including an ordered model for chromosome numbers).
* `visualizeTree()` draws an axis of the pairwise distances used to build
  the tree, as numbers (default) or percentages (`axis = "percent"`), or not
  (`axis = "none"`).  It can mark nodes with points instead of labels
  (`points = TRUE`), and label nothing with `value = NULL`.
* New `makeTraitTibble()` function to copy species traits from a taxon
  table to the tips of a tree, and new `Halo_Taxons` example taxon table.
* The first argument of `visualizeTree()` is renamed from `your_tibble` to
  `tree`.
* New `cladeBars()` function to mark focal clades with bars next to the tips
  of tree plots instead of boxes behind the branches.
* `subTree()` now outputs trees with proper `isTip` and `y` columns.
* Fix `MRCAs()` so that it does not output averages values for species in the
  results table but not in the input tree.
* Fix again accidental discarding of YAML files with `gz` in their name.
* Speed up `formatStats()` roughly 8 times by using `sapply` instead of `do.call`.
* Use `linewidth` instead of the deprecated `label.size` in `visualizeTree()`
  (requires ggplot2 >= 3.5.0).
* Add `hover_text`, `type`, `size` and `alpha` columns to `oikData` `2025_07_17`
  so that it matches the plot data of `MRCA_2D_plot()`.

# ScrambledTreeBuilder 1.3.0

* Expand `recordClades()` to work on taxon tables (one genome per line).
* New `removeAssemblies()` function for on-the-fly filtering.

# ScrambledTreeBuilder 1.2.0

* Fix diagonal values in `Halo_PercentDiff`.
* Fix detection of YAML files.
* Make `formatStats()` robust to the presence of character values.
* Stop dropping zero-valued columns in `formatStats()`.
* New `residualBootstrapTree()` function.
* `makeMatrix()` now warns and returns `NULL` if a column is not found.
* `makeValueTibble()` returns the original tree if `NULL` is passed instead of
   a matrix.
* Added `treeHeatMap()` to visualise distance matrices with their values sorted
  like the branches of a phylogenetic tree, and optional focal clade
  highlighting.

# ScrambledTreeBuilder 1.1.0

* Add `oikData` `2025_02_25` that uses the experimental windowed scrambling
  index.

# ScrambledTreeBuilder 1.0.2

* Make `formatStats()` more robust when computing `percent_aligned`

# ScrambledTreeBuilder 1.0.1

* Fix error messages in `subMatrix()` and document that it also works with
  focal clades.

# ScrambledTreeBuilder 1.0.0

* First release with semantic versioning.

# ScrambledTreeBuilder 0.0.0.9001

* Versions used before June 2025.
