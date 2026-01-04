# Estimate genotyping error rate (REMOVED; will be re-implemented)

Estimate the genotyping error rates in SNP data, based on a pedigree
and/or duplicates. Estimates probabilities (observed given actual)
hom\|other hom, het\|hom, and hom\|het. THESE ARE APPROXIMATE VALUES!

## Usage

``` r
EstEr(
  GenoM,
  Pedigree,
  Duplicates = NULL,
  Er_start = c(0.05, 0.05, 0.05),
  perSNP = FALSE
)
```

## Arguments

- GenoM:

  Genotype matrix

- Pedigree:

  data.frame with columns id - dam - sire

- Duplicates:

  matrix or data.frame with 2 columns, id1 & id2

- Er_start:

  vector of length 3 with starting values for `optim`.

- perSNP:

  logical, estimate error rate per SNP. WARNING not very precise, use
  only as an approximate indicator! Try on simulated data first, e.g.
  with [`SimGeno`](https://jiscah.github.io/reference/SimGeno.md).

## Value

vector of length 3 with estimated genotyping error rates: the
probabilities that

- hom\|hom: an actual homozygote is observed as the other homozygote

- het\|hom: an actual homozygote is observed as heterozygote

- hom\|het: an actual heterozygote is observed as homozygote

These are three independent parameters, that define the genotyping error
matrix (see [`ErrToM`](https://jiscah.github.io/reference/ErrToM.md)) as
follows:

|       |               |            |               |
|-------|---------------|------------|---------------|
|       | **0**         | **1**      | **2**         |
| **0** | \\1-E_1-E_2\\ | \\E_2\\    | \\E_1\\       |
| **1** | \\E_3\\       | \\1-2E_3\\ | \\E_3\\       |
| **2** | \\E_1\\       | \\E_2\\    | \\1-E_1-E_2\\ |

Note that for `optim` a lower bound of 1e-6 and upper bound of 0.499 are
used; if these values are returned this should be interpreted as
'inestimably small' and 'inestimably large', respectively. PLEASE DO NOT
USE THESE VALUES AS INPUT IN SUBSEQUENT ANALYSIS BUT SUBSITUTE BY A
SENSIBLE VALUE!!

## Details

The result should be interpreted as approximate, ballpark estimates! The
estimated error rates from a pedigree will not be as accurate as from
duplicate samples. Errors in individuals without parents or offspring
will not be counted, and errors in individuals with only few offspring
may not be noted either. Deviation of genotype frequencies among
founders from Hardy-Weinberg equilibrium may wrongly be attributed to
genotyping errors. Last but not least, any pedigree errors will result
in higher estimated genotyping errors.

## Examples

``` r
GenoX <- SimGeno(Ped_griffin, nSnp=400, SnpError=c(0.01,0.07, 0.1),
                ParMis=0.1, CallRate=0.9)
# EstEr(GenoM=GenoX, Pedigree=Ped_griffin)
```
