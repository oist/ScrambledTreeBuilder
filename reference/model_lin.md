# Linear interpolation.

This function uses the relationship between cumulative branch length and
time to infer the divergence times of ancestral nodes that have only the
cumulative branch length.

## Usage

``` r
model_lin(known, unknown, x_grid = NULL)
```

## Arguments

- known:

  Nodes with both divergence time (column `divergenceTime`) and
  cumulative branch length (column `cumBranchLen`).

- unknown:

  Nodes with only cumulative branch length (column `cumBranchLen`;
  column `divergenceTime` is `NA`).

- x_grid:

  Optional parameter with extra points beyond those with unknown
  divergence times. These can be used, for example, to plot or check
  accuracy over a wider range.

## Value

A list with 2 types of predictions: 1. `pred_unknown`: predictions for
nodes in the path with unknown divergence time; 2. `pred_grid`:
predictions for the extra points passed in `x_grid`.

## Details

The relationship between cumulative branch length and time is inferred
with Linear interpolation.

## Author

Priscila Biller
