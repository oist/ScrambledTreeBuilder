# Divergence time - TimeTree

Retrieves the divergence time of two species via the TimeTree API.

## Usage

``` r
getDivergenceTime(
  taxId1,
  taxId2,
  timetreeBasePairwise = "http://timetree.temple.edu/api/pairwise/"
)
```

## Arguments

- taxId1:

  An integer: The NCBI id of one taxon

- taxId2:

  An integer: The NCBI id of the other taxon

- timetreeBasePairwise:

  Optional parameter specifying the base URL for retrieving divergence
  times between pairs of species from the TimeTree API.

## Value

A list in which the first item is the divergence time (`divtime`), and
the last two items are the divergence time confidence interval
(`divtimeCI_low` and `divtimeCI_high`). If the divergence time is not
found in the TimeTree database, it returns `NA`.

## Author

Priscila Biller

## Examples

``` r
if (FALSE) { # \dontrun{
getDivergenceTime(9555,9601)} # }
```
