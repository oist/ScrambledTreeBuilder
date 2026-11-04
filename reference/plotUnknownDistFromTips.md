# Plot the relation between distance from tips and lack of information from TimeTree.

Plot the relation between distance from tips and lack of information
from TimeTree.

## Usage

``` r
plotUnknownDistFromTips(unknownTypes, cladeName)
```

## Arguments

- unknownTypes:

  List where each element corresponds to a path from tip to root. Each
  element in the list contains the distances from tips of all nodes
  without TimeTree information for the given path. A distance equals to
  `1` means that the internal node is the direct parent of a leaf,
  whereas a higher distance `n` means that the node is closer to the
  root.

- cladeName:

  The name of the analyzed clade.

## Value

A `ggplot` object.

## Author

Priscila Biller
