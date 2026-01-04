# Simulate Genotyping Errors

Generate errors and missing values in a (simulated) genotype matrix.

## Usage

``` r
MkGenoErrors(
  SGeno,
  CallRate = 0.99,
  SnpError = 5e-04,
  ErrorFV = function(E) c((E/2)^2, E - (E/2)^2, E/2),
  ErrorFM = NULL,
  Error.shape = 0.5,
  CallRate.shape = 1,
  WithLog = FALSE
)
```

## Arguments

- SGeno:

  matrix with genotype data in Sequoia's format: 1 row per individual, 1
  column per SNP, and genotypes coded as 0/1/2.

- CallRate:

  either a single number for the mean call rate (genotyping success), OR
  a vector with the call rate at each SNP, OR a named vector with the
  call rate for each individual. In the third case, ParMis is ignored,
  and individuals in the pedigree (as id or as parent) not included in
  this vector are presumed non-genotyped.

- SnpError:

  either a single value which will be combined with `ErrorFV`, or a
  length 3 vector with probabilities (observed given actual) hom\|other
  hom, het\|hom, and hom\|het; OR a vector or 3XnSnp matrix with the
  genotyping error rate(s) for each SNP.

- ErrorFV:

  function taking the error rate (scalar) as argument and returning a
  length 3 vector with hom-\>other hom, hom-\>het, het-\>hom. May be an
  'ErrFlavour', e.g. 'version2.9'.

- ErrorFM:

  function taking the error rate (scalar) as argument and returning a
  3x3 matrix with probabilities that actual genotype i (rows) is
  observed as genotype j (columns). See below for details. To use, set
  `ErrorFV = NULL`

- Error.shape:

  first shape parameter (alpha) of beta-distribution of per-SNP error
  rates. A higher value results in a flatter distribution.

- CallRate.shape:

  as Error.shape, for per-SNP call rates.

- WithLog:

  Include dataframe in output with which datapoints have been edited,
  with columns id - SNP - actual (original, input) - observed (edited,
  output).

## Value

The input genotype matrix, with some genotypes replaced, and some set to
missing (-9). If `WithLog=TRUE`, a list with 3 elements: GenoM, Log, and
Counts_actual (genotype counts in input, to allow double checking of
simulated genotyping error rate).
