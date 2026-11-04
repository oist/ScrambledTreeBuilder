# Cumulative branch lengths

This method computes the cumulative sum of branch lengths along the path
from the root to each node in the tree.

## Usage

``` r
computeCumulativeBranchLengths(tree, rootId = NULL)
```

## Arguments

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

- rootId:

  The node ID of the root (column `node` in `tree`).

## Value

The given tree with an additional `cumBranchLen` column containing the
cumulative sum of branch lengths for each node.

## Author

Priscila Biller
