# Back-transform IDs

Reverse the joining of FID and IID in
[`GenoConvert`](https://jiscah.github.io/reference/GenoConvert.md) and
[`LHConvert`](https://jiscah.github.io/reference/LHConvert.md)

## Usage

``` r
PedStripFID(Ped, FIDsep = "__")
```

## Arguments

- Ped:

  pedigree as returned by sequoia (e.g. `SeqOUT$Pedigree`).

- FIDsep:

  characters inbetween FID and IID in composite-ID.

## Value

A pedigree with 6 columns

- FID:

  family ID of focal individual (offspring).

- id:

  within-family of focal individual

- dam.FID:

  original family ID of assigned dam

- dam:

  within-family of dam

- sire.FID:

  original family ID of assigned sire

- sire:

  within-family of sire

## Details

Note that the family IDs are the ones provided, and not automatically
updated. New, numeric ones can be obtained with
[`FindFamilies`](https://jiscah.github.io/reference/FindFamilies.md).
