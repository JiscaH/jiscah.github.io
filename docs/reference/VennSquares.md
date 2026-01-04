# Square Venn diagram

Draw Venn diagram with squares, with match/mismatch in overlapping area.

## Usage

``` r
VennSquares(count, BL = c(0, 0), COL, withText = TRUE, withLegend = FALSE)
```

## Arguments

- count:

  a length 5 named vector: 'Total', 'Match', 'Mismatch', 'P1only', and
  'P2only'.

- BL:

  a length 2 vector with coordinates of bottom-mid of Ped1 square.

- COL:

  a length 4 character vector with colours, named 'Match', 'Mismatch',
  'Ped1', 'Ped2'.

- withText:

  logical, add count to each rectangle.

- withLegend:

  logical, add legend at the bottom of the plot.

## See also

[`PlotPedComp`](https://jiscah.github.io/reference/PlotPedComp.md)
