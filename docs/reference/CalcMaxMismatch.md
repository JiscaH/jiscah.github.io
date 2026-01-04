# Maximum Number of Mismatches

Calculate the maximum expected number of mismatches for duplicate
samples, parent-offspring pairs, and parent-parent-offspring trios.

## Usage

``` r
CalcMaxMismatch(
  Err,
  MAF,
  ErrFlavour = "version2.9",
  qntl = 1 - 1e-05,
  Return = "Counts"
)
```

## Arguments

- Err:

  estimated genotyping error rate, as a single number or 3x3 matrix
  (averaged value(s) across SNPs), or a vector with the same length as
  MAF, or a nSnp x 3 x 3 array. If a matrix, this should be the
  probability of observed genotype (columns) conditional on actual
  genotype (rows). Each row must therefore sum to 1. If an array, each
  3x3 slice should abide this rule.

- MAF:

  vector with minor allele frequency at each SNP.

- ErrFlavour:

  function that takes `Err` as input, and returns a 3x3 matrix of
  observed (columns) conditional on actual (rows) genotypes, or choose
  from inbuilt ones as used in sequoia 'version2.0', 'version1.3', or
  'version1.1'. Ignored if `Err` is a matrix. See
  [`ErrToM`](https://jiscah.github.io/reference/ErrToM.md).

- qntl:

  quantile of binomial distribution to be used as the maximum, of
  individual-level probability. For a desired dataset-level probability
  quantile \\Q\\, use `qntl`\\= Q^{(1/N)}\\, where \\N\\ is the number
  of individuals.

- Return:

  Either 'Counts' to return the threshold counts (default), or 'Probs'
  to return the mismatch probabilities from which these counts are
  calculated.

## Value

A vector with three integers:

- DUP:

  Maximum number of differences between 2 samples from the same
  individual

- OH:

  Maximum number of Opposing Homozygous SNPs between a true
  parent-offspring pair

- ME:

  Maximum number of Mendelian Errors among a true parent-parent-
  offspring trio

.

## Details

The thresholds for maximum number of mismatches calculated here aim to
minimise false negatives, i.e. to minimise the chance that any true
duplicates or true parent-offspring pairs are already excluded during
the filtering steps where these `MaxMismatch` values are used.
Consequently, there is a high probability of false positives, i.e. it is
likely that some sample pairs with fewer mismatches than the
`MaxMismatch` threshold, are in fact not duplicate samples or
parent-offspring pairs. Use of these `MaxMismatch` thresholds is
therefore only the first step of pedigree reconstruction by
[`sequoia`](https://jiscah.github.io/reference/sequoia.md).

## See also

[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md).

## Examples

``` r
CalcMaxMismatch(Err = 0.05, MAF = runif(n=100, min=0.3, max=0.5))
#> DUP  OH  ME 
#>  24   8  14 

# in sequoia() qntl depends on the number of genotyped individuals, to get an
# approximately constant false exclusion rate at dataset-level
sts <- SnpStats(Geno_griffin, Plot=FALSE, quiet=TRUE, calc_HWE=FALSE)
MAF <- ifelse(sts[,'AF'] < 0.5, sts[,'AF'], 1-sts[,'AF'])
sequoia::CalcMaxMismatch(Err = 0.001,
                         MAF = MAF,
                         qntl = 0.9999^(1/nrow(Geno_griffin)))
#> DUP  OH  ME 
#>   8   4   6 
```
