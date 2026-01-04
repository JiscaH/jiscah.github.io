# LLR-age from Ageprior Matrix

Get log10-likelihood ratios for a specific age difference from matrix
`AgePriorExtra`.

## Usage

``` r
GetLLRAge(AgePriorExtra, agedif, patmat)
```

## Arguments

- AgePriorExtra:

  matrix in [`sequoia`](https://jiscah.github.io/reference/sequoia.md)
  output

- agedif:

  vector with age differences, in whole numbers. Must occur in rownames
  of `AgePriorExtra`.

- patmat:

  numeric vector; choose maternal (1), paternal (2) relatives, or for
  each relationship the most-likely alternative (3).

## Value

A matrix with `nrow` equal to the length of `agedif`, and 7 columns:
PO-FS-HS-GP-FA-HA-U.

## Details

This is a simple helper function to extract values from `AgePriorExtra`,
e.g. to use together with
[`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md).

## Examples

``` r
# For a pair with unknown age difference, explore the difference age-based
# LLRs for all relationships, for a range of plausible age differences.
PairsG <- data.frame(ID1 = 'A', ID2 = 'B', AgeDif = rep(c(-2,2,3),2),
                   PatMat = rep(1:2, each=3))
cbind(PairsG,
      GetLLRAge(SeqOUT_griffin$AgePriorExtra,
                agedif = PairsG$AgeDif, patmat = PairsG$PatMat))
#>   ID1 ID2 AgeDif PatMat    PO    FS    HS    GP    FA    HA U
#> 1   A   B     -2      1  -Inf -0.65 -0.80  -Inf -1.49 -1.20 0
#> 2   A   B      2      1  0.42 -0.65 -0.80  0.14  0.25  0.28 0
#> 3   A   B      3      1 -0.96 -0.95 -1.10  0.44  0.06 -0.06 0
#> 4   A   B     -2      2  -Inf -0.65 -0.13  -Inf -1.72 -1.43 0
#> 5   A   B      2      2  0.52 -0.65 -0.13 -0.09  0.21  0.28 0
#> 6   A   B      3      2 -0.09 -0.95 -0.43  0.35  0.18  0.10 0
```
