# Fortran Simulate Genotyping Errors

Wrapper for Fortran function to simulate genotyping errors.

## Usage

``` r
DoErrors(SGeno, Act2Obs)
```

## Arguments

- SGeno:

  matrix with genotype data, size nInd x nSnp.

- Act2Obs:

  array with conditional probability of observing genotype i conditional
  on actual genotype j, size nSnp x 3 x 3.

## Value

`SGeno` with errors.
