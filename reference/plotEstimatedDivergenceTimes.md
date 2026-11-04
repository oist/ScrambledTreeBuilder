# Plot the divergence times estimated with different models by extrapolation.

Plot the divergence times estimated with different models by
extrapolation.

## Usage

``` r
plotEstimatedDivergenceTimes(cladeName, estimatedResults, expectedVals = NULL)
```

## Arguments

- cladeName:

  The name of the analyzed clade.

- estimatedResults:

  The output from the function `inferDivTimesFromExtrapolation`.

- expectedVals:

  An optional dataframe with columns `cumBranchLen` and `divergenceTime`
  for checking purposes.

## Value

A list with the following: 1. `plot` : A `ggplot` object. 2. `plotLabel`
: A title for the plot.

## Author

Priscila Biller
