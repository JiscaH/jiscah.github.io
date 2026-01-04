# Dummifiable IDs

Get the dummifiable individuals, using various possible criteria

## Usage

``` r
GetDummifiable(Pedigree, gID, minSibSize)
```

## Arguments

- Pedigree:

  dataframe with id - dam - sire.

- gID:

  vector with IDs of SNP-genotyped individuals.

- minSibSize:

  minimum requirements to be considered dummifiable:

  - '1sib' : sibship of size 1, with or without grandparents. The latter
    aren't really a sibship, but can be useful in some situations.

  - '1sib1GP': sibship of size 1 with at least 1 grandparent

  - '2sib': at least 2 siblings, with or without grandparents. Used by
    [`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)

  .

## Value

A length-2 list (dams, sires) with each element a vector with
dummifiable ids

## Details

values of minSibSize used by calling functions

- 1sib:

  CalcOHLLR, CalcPairLL

- 1sib1GP:

  getAssignCat (default when user called)

- 2sib:

  PedCompare
