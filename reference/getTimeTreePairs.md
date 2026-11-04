# Gets all pairs of nodes with TimeTree info for all nodes.

Compute all pairs of leaves that could provide information about the
divergence time of an ancestral node.

## Usage

``` r
getTimeTreePairs(nodes, tree, ncbi_ids_timetree)
```

## Arguments

- nodes:

  A list in which each element corresponds to an ancestral node. Each
  element was computed with the function `getLeaves` and consists of a
  named list that maps the ID of each direct descendant node to the IDs
  of its descendant leaves.

- tree:

  A
  [`ConvenientTblTree`](https://brennern.github.io/ScrambledTreeBuilder/reference/ConvenientTblTree.md)
  object.

- ncbi_ids_timetree:

  A list of NCBI ids that appear in the TimeTree database.

## Value

A named list that maps each ancestral node to pairs of leaves,
identified by their node IDs. These pairs of leaves are informative for
estimating the divergence time of that ancestor.

## Author

Priscila Biller
