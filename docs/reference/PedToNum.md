# Turn Character Pedigree into Numeric Pedigree

Genotyped individuals get rownumber in genotype matrix, non-genotyped
individuals either all get an arbitrary negative number
(`DoDummies = 'new'`) or only individuals with a dummy ID get the
corresponding negative number (`DoDummies = 'old'`). Note that the
number series will overlap for dummy males and dummy females.

## Usage

``` r
PedToNum(
  Pedigree = NULL,
  gID = NULL,
  DoDummies = "new",
  DumPrefix = c("F0", "M0")
)
```

## Arguments

- Pedigree:

  dataframe with id - dam - sire. It is assumed
  [`PedPolish`](https://jiscah.github.io/reference/PedPolish.md) has
  been called beforehand so that column names are correct and all
  columns are as.character.

- gID:

  vector with IDs of SNP-genotyped individuals.

- DoDummies:

  'new', 'old', or 'no' (ignore all non-genotyped individuals).

- DumPrefix:

  Prefix to identify dummies when `DoDummies = 'old'`

## Value

a list with

- PedPar:

  An nInd x 2 matrix with the numeric IDs of parents of genotyped
  individuals

- DumPar:

  A matrix with parents of dummies, see
  [`FoldSibGPs`](https://jiscah.github.io/reference/FoldSibGPs.md)

- Renamed:

  a length-2 list (dams, sires) with each element a dataframe with
  columns: 'name' (original character ID), 'num' (number ID, negative)
  for each dummified individual

- Nd:

  a length 2 vector, no. dummies found/created for dams and sires

## Details

If `DoDummies='new'`,
[`getAssignCat`](https://jiscah.github.io/reference/getAssignCat.md) is
used with `minSibSize ="1sib1GP"`, and any existing dummy coding is
ignored (F0001, F0002 may become -3, -6). If `DoDummies='old'`, the
existing dummy coding is respected (F0001, F0002 will become -1, -2),
but other non-genotyped individuals are ignored.
