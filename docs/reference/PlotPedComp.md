# Visualise PedCompare Output

square Venn diagrams with
[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)
`Counts`.

## Usage

``` r
PlotPedComp(Counts, sameSize = FALSE)
```

## Arguments

- Counts:

  a 7x5x2 array with counts of matches and mismatches per category
  (genotyped vs dummy), as returned by
  [`PedCompare`](https://jiscah.github.io/reference/PedCompare.md).

- sameSize:

  logical, make all per-category Venn diagrams the same size `TRUE`, or
  make their size proportional to the counts (`FALSE`, the default). If
  `TRUE`, a warning is printed at the bottom.

## See also

[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)

## Examples

``` r
PC.g <- PedCompare(Ped1 = cbind(FieldMums_griffin, sire=NA),
                   Ped2 = SeqOUT_griffin$Pedigree)

PlotPedComp(PC.g$Counts)
```
