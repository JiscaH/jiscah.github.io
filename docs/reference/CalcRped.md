# Calculate Pedigree Relatedness

Morph pedigree into a kinship2 compatible format and use
[`kinship`](https://rdrr.io/pkg/kinship2/man/kinship.html) to calculate
kinship coefficients; relatedness = 2\*kinship.

## Usage

``` r
CalcRped(Pedigree, OUT = "DF")
```

## Arguments

- Pedigree:

  dataframe with columns id-dam-sire.

- OUT:

  desired output format, 'M' for matrix or 'DF' for dataframe with
  columns IID1 - IID2 - R.ped.

## Value

A matrix or dataframe.
