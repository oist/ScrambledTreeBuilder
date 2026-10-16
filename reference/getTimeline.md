# Timeline - TimeTree

Fetches a list of divergence times of nodes in TimeTree from the
specified taxon to the common ancestor of all cellular organisms.

## Usage

``` r
getTimeline(
  taxId,
  divergenceTimeType = "adjusted_age",
  timetreeBaseTimeline = "http://timetree.temple.edu/api/timeline/"
)
```

## Arguments

- taxId:

  An integer: The NCBI id of a taxon

- divergenceTimeType:

  Optional parameter specifying which type of divergence time to read
  from TimeTree. Two types are available; the default is `adjusted_age`:

  1.  `summary_age`: the median divergence time across all studies from
      TimeTree;

  2.  `adjusted_age`: a time adjusted when a node is older than its
      parent in the global timetree reconstructed from individual
      timetrees.

- timetreeBaseTimeline:

  Optional parameter specifying the base URL for retrieving timelines
  from the TimeTree API.

## Value

A dataframe in which each row corresponds to an ancestor. Each row
contains the divergence time (`divtime`), its confidence interval
(`divtimeCI_low` and `divtimeCI_high`), and additional information,
including the scientific name (`name`) and taxonomic rank (`rank`). If
the divergence time is not found in the TimeTree database, it returns
`NA`.

## Author

Priscila Biller

## Examples

``` r
if (FALSE) { # \dontrun{
getTimeline(9555)} # }
```
