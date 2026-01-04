# Estimate Genotyping Error Rate

Estimate genotyping error rate from Mendelian errors per SNP.

## Usage

``` r
OHperSNP(GenoM, Par, Dups = NULL)
```

## Arguments

- GenoM:

  genotype matrix, in sequoia's format: 1 column per SNP, 1 row per
  individual, genotypes coded as 0/1/2/-9, and rownames giving
  individual IDs.

- Par:

  pedigree dataframe, only genotyped parents are used.

- Dups:

  pairs of duplicates

## Value

A dataframe with columns:

- n.dam, n.sire, n.pair:

  Number of dams, sires, parent-pairs successfully genotyped for the SNP

- OHdam, OHsire:

  Count of number of opposing homozygous cases

- MEpair:

  Count of Mendelian errors, includes opposing homozygous cases

- n.dups, n.diff:

  Number of duplicate pairs successfully genotyped for the SNP; number
  of differences

## See also

[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md).
