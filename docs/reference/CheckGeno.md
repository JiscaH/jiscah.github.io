# Check Genotype Matrix

Check that the provided genotype matrix is in the correct format, and
check for low call rate samples and SNPs.

## Usage

``` r
CheckGeno(
  GenoM,
  quiet = FALSE,
  Plot = FALSE,
  Return = "GenoM",
  Strict = TRUE,
  DumPrefix = c("F0", "M0")
)
```

## Arguments

- GenoM:

  the genotype matrix.

- quiet:

  suppress messages.

- Plot:

  display the plots of
  [`SnpStats`](https://jiscah.github.io/reference/SnpStats.md).

- Return:

  either 'GenoM' to return the cleaned-up genotype matrix, or 'excl' to
  return a list with excluded SNPs and individuals (see Value).

- Strict:

  Exclude any individuals genotyped for \<5 genotyped for \<5 up to
  version 2.4.1. Otherwise only excluded are (very nearly) monomorphic
  SNPs, SNPs scored for fewer than 2 individuals, and individuals scored
  for fewer than 2 SNPs.

- DumPrefix:

  length 2 vector, to check if these don't occur among genotyped
  individuals.

## Value

If `Return='excl'` a list with, if any are found:

- ExcludedSNPs:

  SNPs scored for \<10 excluded when running
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md)

- ExcludedSnps-mono:

  monomorphic (fixed) SNPs; automatically excluded when running
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md). This
  includes nearly-fixed SNPs with MAF \\= 1/2N\\. Column numbers are
  \*after\* removal of `ExcludedSNPs`, if any.

- ExcludedIndiv:

  Individuals scored for \<5 reliably included during pedigree
  reconstruction. Individual call rate is calculated after removal of
  'Excluded SNPs'

- Snps-LowCallRate:

  SNPs scored for 10 recommended to be filtered out

- Indiv-LowCallRate:

  individuals scored for \<50 recommended to be filtered out

When `Return='excl'` the return is
[`invisible`](https://rdrr.io/r/base/invisible.html), i.e. a check is
run and warnings or errors are always displayed, but nothing may be
returned.

## Thresholds

Appropriate call rate thresholds for SNPs and individuals depend on the
total number of SNPs, distribution of call rates, genotyping errors, and
the proportion of candidate parents that are SNPd (sibship clustering is
more prone to false positives). Note that filtering first on SNP call
rate tends to keep more individuals in.

## See also

[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md) to
calculate SNP call rates;
[`CalcOHLLR`](https://jiscah.github.io/reference/CalcOHLLR.md) to count
the number of SNPs scored in both focal individual and parent.

## Examples

``` r
GenoM <- SimGeno(Ped_HSg5, nSnp=400, CallRate = runif(400, 0.2, 0.8))
# the quick way:
GenoM.checked <- CheckGeno(GenoM, Return="GenoM")
#> ! There are 191 SNPs scored for <50% of individuals
#> ℹ There are  1000  individuals and  400  SNPs.

# the user supervised way:
Excl <- CheckGeno(GenoM, Return = "excl")
#> ! There are 191 SNPs scored for <50% of individuals
#> ℹ There are  1000  individuals and  400  SNPs.
GenoM.orig <- GenoM   # make a 'backup' copy
if ("ExcludedSnps" %in% names(Excl))
  GenoM <- GenoM[, -Excl[["ExcludedSnps"]]]
if ("ExcludedSnps-mono" %in% names(Excl))
  GenoM <- GenoM[, -Excl[["ExcludedSnps-mono"]]]
if ("ExcludedIndiv" %in% names(Excl))
  GenoM <- GenoM[!rownames(GenoM) %in% Excl[["ExcludedIndiv"]], ]

# warning about  SNPs scored for <50% of individuals ?
# note: this is not necessarily a problem, and sometimes unavoidable.
SnpCallRate <- apply(GenoM, MARGIN=2,
                     FUN = function(x) sum(x!=-9)) / nrow(GenoM)
hist(SnpCallRate, breaks=50, col="grey")

GenoM <- GenoM[, SnpCallRate > 0.6]

# to filter out low call rate individuals: (also not necessarily a problem)
IndivCallRate <- apply(GenoM, MARGIN=1,
                       FUN = function(x) sum(x!=-9)) / ncol(GenoM)
hist(IndivCallRate, breaks=50, col="grey")

GoodSamples <- rownames(GenoM)[ IndivCallRate > 0.8]
```
