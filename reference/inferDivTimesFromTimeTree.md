# Divergence time based on TimeTree

Infer the divergence time of one node based on the TimeTree information.

## Usage

``` r
inferDivTimesFromTimeTree(node_pairs, timelines)
```

## Arguments

- node_pairs:

  A matrix with two columns (`leaf_1` and `leaf_2`) computed with the
  function `getTimeTreePairs`. Each row in the matrix corresponds to a
  pair of informative leaves, identified by their node ID.

- timelines:

  A named list mapping each node to its timeline from TimeTree. Each
  timeline contains all ancestors of the node up to the root of the Tree
  of Life. A timeline is represented as a dataframe with one row per
  ancestor.

## Value

A data frame with one row per informative pair of leaves. Each row is
derived from the `timelines` object and contains the node IDs of the
paired leaves, along with information about their most recent common
ancestor (MRCA), including its divergence time.

## Details

This method estimates the divergence time of a node given "informative"
pairs of descendant leaves.

An informative pair of leaves consists of two leaves, each belonging to
a different direct descendant of the node, for which TimeTree data are
available.

For each pair, the method identifies the most recent common ancestor
(MRCA) and retrieves its estimated divergence time. Ideally, all pairs
should have the same MRCA (i.e., the given node) and yield the same
divergence time estimate. If they do not, this indicates that the
topology of the phylogeny from TimeTree differs from that inferred by
the ScrambledTreeBuilder package.

## Author

Priscila Biller
