# Plot relation between divergence time and cumulative branch length for a clade.

This function creates a plot that shows the relationship between
divergence time (taken from TimeTree) and cumulative branch length for
the whole clade.

## Usage

``` r
plotDivergenceTimesFromTimeTree(tree, cladeName)
```

## Arguments

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

- cladeName:

  The name of the analyzed clade.

## Value

A `ggplot` object.

## Author

Priscila Biller
