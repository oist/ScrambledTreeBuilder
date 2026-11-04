# Gets the most recent common ancestor (MRCA) of two nodes.

Computes the most recent common ancestor of two nodes, given their
timelines from TimeTree.

## Usage

``` r
getMRCA(taxid1, taxid2, timelines)
```

## Arguments

- taxid1:

  An integer: The NCBI id of one taxon

- taxid2:

  An integer: The NCBI id of the other taxon

- timelines:

  A named list mapping each node to its timeline from TimeTree. Each
  timeline contains all ancestors of the node up to the root of the Tree
  of Life. A timeline is represented as a dataframe with one row per
  ancestor.

## Value

The row in the dataframe corresponding to the most recent common
ancestor of the two given nodes.

## Author

Priscila Biller
