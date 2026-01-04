# Make Pairs Fortran Compatible

Convert dataframe `Pairs` into a list of integer vectors. Called only by
[`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md).

## Usage

``` r
FortifyPairs(Pairs, gID, Renamed, LH, Ped)
```

## Arguments

- Pairs:

  dataframe with columns ID1 - ID2 - Sex1 - Sex2 - AgeDif - focal - k.

- gID:

  character vector with IDs of genotyped individuals.

- Renamed:

  length-2 list (dams, sires) each with a 2-column dataframe. matching
  character IDs to negative numbers, for dummified individuals. Element
  of the list returned by
  [`PedToNum`](https://jiscah.github.io/reference/PedToNum.md).

- LH:

  lifehistory dataframe, ID - Sex - BirthYear.

- Ped:

  pedigree, to ensure dams have sex=1 & sires sex=2

## Value

A named list, with elements ID - Sex - AgeDif - focal. The first two are
per individual and thus each have length 2\*nrow(Pairs), while the last
two have length 1\*nrow(Pairs).
