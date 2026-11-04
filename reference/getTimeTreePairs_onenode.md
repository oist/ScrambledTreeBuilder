# Gets all pairs of nodes with TimeTree info for a single node.

Computes all pairs of leaves that could provide information about the
divergence time of one ancestral node.

## Usage

``` r
getTimeTreePairs_onenode(children)
```

## Arguments

- children:

  A named list of all (informative) leaves descending from a given node,
  grouped by its direct descendants. The names of the elements are the
  node IDs of the direct descendants.

## Value

A matrix in which each row corresponds to a pair of leaves.

## Author

Priscila Biller
