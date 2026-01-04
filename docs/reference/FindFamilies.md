# Assign Family IDs

Find clusters of connected individuals in a pedigree, and assign each
cluster a unique family ID (FID).

## Usage

``` r
FindFamilies(Pedigree = NULL, SeqList = NULL, MaybeRel = NULL)
```

## Arguments

- Pedigree:

  dataframe with columns id - parent1 - parent2; only the first 3
  columns will be used.

- SeqList:

  list as returned by
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md). If
  `Pedigree` is not provided, the element `Pedigree` from this list will
  be used if present, and element `Pedigreepar` otherwise.

- MaybeRel:

  Output from
  [`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md), a
  dataframe with probable but non-assigned relatives.

## Value

A numeric vector with length equal to the number of unique individuals
in the pedigree (i.e. number of rows in pedigree after running
[`PedPolish`](https://jiscah.github.io/reference/PedPolish.md) on
`Pedigree`).

## Details

This function repeatedly finds all ancestors and all descendants of each
individual in turn, and ensures they all have the same Family ID. Not
all connected individuals are related, e.g. all grandparents of an
individual will have the same FID, but will typically be unrelated.

When `UseMaybeRel = TRUE`, probable relatives are added to existing
family clusters, or existing family clusters may be linked together.
Currently no additional family clusters are created.

## See also

[`GetAncestors`](https://jiscah.github.io/reference/GetAncestors.md)`, `[`GetDescendants`](https://jiscah.github.io/reference/GetDescendants.md)`, `[`getGenerations`](https://jiscah.github.io/reference/getGenerations.md)

## Examples

``` r
PedG <- SeqOUT_griffin$PedigreePar[,1:3]
FID_G <- FindFamilies(PedG)
PedG[FID_G==4,]
#>             id         dam        sire
#> 6  i014_2001_F        <NA>        <NA>
#> 10 i018_2001_M        <NA>        <NA>
#> 15 i025_2002_M i014_2001_F i018_2001_M
#> 37 i057_2003_M        <NA> i018_2001_M
#> 53 i078_2004_M        <NA> i025_2002_M
#> 64 i094_2005_M        <NA> i078_2004_M
#> 93 i135_2007_F        <NA> i078_2004_M
```
