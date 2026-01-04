# Chi-square Test on Observed vs Expected Genotypes

For one SNP and all offspring-parent-parent trios or single
parent-offspring pairs, calculate the expected genotype frequencies
given the allele frequency, genotyping error rate, and error flavour,
and perform a chi-square test.

## Usage

``` r
CalcChi2(E, q, A.obs, ErrF)
```

## Arguments

- E:

  presumed genotyping error rate.

- q:

  allele frequency.

- A.obs:

  array of dim 4x4x4 with counts of joined offspring-parent-parent at
  the SNPs

- ErrF:

  ErrFlavour; function that takes the genotyping error rate `Err` as
  input, and returns a 3x3 matrix of observed (columns) conditional on
  actual (rows) genotypes, or choose from inbuilt ones as used in
  sequoia 'version2.0', 'version1.3', or 'version1.1'. See
  [`ErrToM`](https://jiscah.github.io/reference/ErrToM.md).

## Value

The chisquare value of the test.

## See also

[`EstErr`](https://jiscah.github.io/reference/EstErr.md)`, `[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md).

## Examples

``` r
if (FALSE) {
E.hat <- numeric(ncol(GenoM))  # based on trio where possible
for (l in 1:ncol(GenoM)) {
    E.hat[l] <- stats::optimise(CalcChi2, interval=c(0,1), q=AF[l],
                         A.obs=OO.trio[l,,,], ErrF=ErrF)$minimum
} }
```
