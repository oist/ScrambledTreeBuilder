# Average results

Most of the times we have two sets of alignment statistics for a pair of
genomes, where the *target* and *query* roles are exchanged. Sometimes
we want to plot only one point per pair. This function produces the
averages.

## Usage

``` r
averageResults(df)
```

## Arguments

- df:

  A results data frame produced by
  [`formatStats()`](https://brennern.github.io/ScrambledTreeBuilder/reference/formatStats.md)
  (therefore) having a proper `lab` column.

## Value

A tibble, even when the input is a data frame, because we use *dplyr*…

## Author

Charles Plessy

## Examples

``` r
averageResults(Halo_DF)
#> Warning: There was 1 warning in `summarise()`.
#> ℹ In argument: `across(where(is.numeric), mean, na.rm = TRUE)`.
#> ℹ In group 1: `lab = "Halobacterium_noricense\nHalobacterium_litoreum"`.
#> Caused by warning:
#> ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
#> Supply arguments directly to `.fns` through an anonymous function instead.
#> 
#>   # Previously
#>   across(a:b, mean, na.rm = TRUE)
#> 
#>   # Now
#>   across(a:b, \(x) mean(x, na.rm = TRUE))
#> ℹ The deprecated feature was likely used in the ScrambledTreeBuilder package.
#>   Please report the issue at
#>   <https://github.com/brennern/ScrambledTreeBuilder/issues>.
#> # A tibble: 15 × 237
#>    lab                aligned_length_Min aligned_length_Q1 aligned_length_Median
#>    <chr>                           <dbl>             <dbl>                 <dbl>
#>  1 "Halobacterium_no…               64                718.                 1416.
#>  2 "Halobacterium_sa…               87                778.                 1457.
#>  3 "Halobacterium_sa…               76                724.                 1345.
#>  4 "Haloferax_medite…               84                562.                  972.
#>  5 "Haloferax_medite…               76.5              617                  1042 
#>  6 "Haloferax_medite…               77                591.                  964.
#>  7 "Haloferax_volcan…               72.5              586.                  996 
#>  8 "Haloferax_volcan…               79                586                  1001.
#>  9 "Haloferax_volcan…               77                626                  1001.
#> 10 "Haloferax_volcan…               74                733.                 1580.
#> 11 "Salarchaeum_japo…               72.5              673.                 1256.
#> 12 "Salarchaeum_japo…               77                680.                 1236.
#> 13 "Salarchaeum_japo…               79                631                  1104 
#> 14 "Salarchaeum_japo…               73                566.                  942.
#> 15 "Salarchaeum_japo…               76                622.                 1004.
#> # ℹ 233 more variables: aligned_length_Mean <dbl>, aligned_length_Q3 <dbl>,
#> #   aligned_length_Max <dbl>, aligned_length_L50 <dbl>,
#> #   aligned_length_Total <dbl>, aligned_length_N <dbl>,
#> #   aligned_score_Min <dbl>, aligned_score_Q1 <dbl>,
#> #   aligned_score_Median <dbl>, aligned_score_Mean <dbl>,
#> #   aligned_score_Q3 <dbl>, aligned_score_Max <dbl>, aligned_score_L50 <dbl>,
#> #   aligned_score_Total <dbl>, aligned_score_N <dbl>, …
```
