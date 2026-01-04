# wrapper for Fortran function to calculate total likelihood

function to be optimised

## Usage

``` r
CalcLL(Er_hat, GenoM, Parents, DupsV)
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

the negative of the total log10-likelihood
