# Corner coordinates

Calculate corner coordinates for each of the four rectangles in a square
Venn diagram

## Usage

``` r
CalcCorners(count)
```

## Arguments

- count:

  a length 5 named vector: 'Total', 'Match', 'Mismatch', 'P1only', and
  'P2only'.

## Value

a 4x4 matrix with columns "xleft", "xright", "ybottom", "ytop" (as used
by [`rect`](https://rdrr.io/r/graphics/rect.html)) and rows "Ped1",
"Ped2", "Mismatch", "Match".

## Details

the bottom-left corner of the Ped1 square is (0,0); offset is done by
[`VennSquares`](https://jiscah.github.io/reference/VennSquares.md). The
size of the Ped1 and Ped2 squares is proportional to their count, i.e.
N1 = count\["Total"\] - count\["P2only"\], and the length of each size
thus proportional to the `sqrt` of that.

The x-location of the Ped2 square is a function of the amount of overlap
(Match + Mismatch): if 0 coco\["Ped1", "xright"\], if 100 coco\["Ped1",
"xright"\]; and proportional in-between these two extremes.

The overlap area between Ped1 and Ped2 is split into Mismatch (bottom)
and Match (top).

## See also

[`PlotPedComp`](https://jiscah.github.io/reference/PlotPedComp.md)`, `[`VennSquares`](https://jiscah.github.io/reference/VennSquares.md)
