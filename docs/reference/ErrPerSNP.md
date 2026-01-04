# wrapper for Fortran function to estimate error rate for each SNP

WARNING: not very precise, especially not for low error rates.

## Usage

``` r
ErrPerSNP(Er_hat, GenoM, Parents, DupsV)
```

## Arguments

- Er_hat:

  length 3 vector with genotyping error rates

- GenoM:

  Genotype matrix

- Parents:

  Pedigree, already converted to rownumbers in GenoM

- DupsV:

  vector with duplicate samples, already converted to rownumbers

## Value

a matrix with 3 columns: for each SNP the probabilities (observed given
actual) hom\|other hom, het\|hom, and hom\|het.
