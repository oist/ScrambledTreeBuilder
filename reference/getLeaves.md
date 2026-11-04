# Leaves grouped by direct descendants of a given node

Get all leaves descending from a given node, grouped by its direct
descendants.

## Usage

``` r
getLeaves(node_id, tree)
```

## Arguments

- node_id:

  An integer corresponding to the node ID.

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

## Value

A named list mapping each direct descendant node ID to the IDs of its
descendant leaves (a list of lists).

## Author

Priscila Biller
