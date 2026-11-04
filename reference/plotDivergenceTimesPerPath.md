# Plot stats about each path from tip to root.

This function creates a plot that shows, for each path from tip to root:

1.  how many ancestral nodes have information from TimeTree;

2.  how many ancestral nodes do not have information from TimeTree
    because all tips belong to the same species;

3.  how many ancestral nodes do not have information from TimeTree, but
    tips are different species.

## Usage

``` r
plotDivergenceTimesPerPath(tree, pathEstimates, cladeName)
```

## Arguments

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

- pathEstimates:

  Named list of dataframes, where each dataframe corresponds to a path
  from tip to root and its results.

- cladeName:

  The name of the analyzed clade.

## Value

A `ggplot` object.

## Author

Priscila Biller
