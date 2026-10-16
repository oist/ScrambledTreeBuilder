# Convenient tbl_trees

Wrapper S7 class for the `tbl_tree` class of the `tidytree` package.

## Usage

``` r
ConvenientTblTree
```

## Details

*Convenient tbl_tree* objects are typically constructed with the
[`makeTidyTree()`](https://brennern.github.io/ScrambledTreeBuilder/reference/makeTidyTree.md)
function. In addition to usual `tbl_tree` objects, they also have a
`isTip` column (boolean) and a `y` column that indicates the order in
which the nodes would be plotted. This is needed by the
[`orderWithTree()`](https://brennern.github.io/ScrambledTreeBuilder/reference/orderWithTree.md)
function that helps making heatmaps sorted in the same way.

## Author

Charles Plessy
