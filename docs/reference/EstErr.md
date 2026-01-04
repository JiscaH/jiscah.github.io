# Estimate Genotyping Error Rate

Estimate genotyping error rate from Mendelian errors per SNP.

## Usage

``` r
EstErr(GenoM, Par, ErrFlavour = "version2.0")
```

## Arguments

- GenoM:

  genotype matrix, in sequoia's format: 1 column per SNP, 1 row per
  individual, genotypes coded as 0/1/2/-9, and rownames giving
  individual IDs.

- Par:

  pedigree dataframe, only genotyped parents are used.

- ErrFlavour:

  function that takes the genotyping error rate `Err` as input, and
  returns a 3x3 matrix of observed (columns) conditional on actual
  (rows) genotypes, or choose from inbuilt ones as used in sequoia
  'version2.0', 'version1.3', or 'version1.1'. See
  [`ErrToM`](https://jiscah.github.io/reference/ErrToM.md).

## Value

A dataframe with columns:

- Err.hat:

  Error rate, as estimated from the joined offspring-parent (-parent)
  genotypes and the presumed error structure (`ErrFlavour`)

- n.dam, n.sire, n.pair:

  Number of dams, sires, parent-pairs succesfully genotyped for the SNP

- OHdam, OHsire:

  Count of number of opposing homozygous cases

- MEpair:

  Count of Mendelian errors, includes opposing homozygous cases

## See also

[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md).
