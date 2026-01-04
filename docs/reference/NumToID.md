# Change Numeric Pedigree back to Character Pedigree

Reverse [`PedToNum`](https://jiscah.github.io/reference/PedToNum.md), 1
column at a time.

## Usage

``` r
NumToID(x, k = 0, gID = NULL, DumPrefix = c("F", "M"))
```

## Arguments

- x:

  vector with numbers.

- k:

  1=dam, 2=sire, needed to distinguish dummy females from dummy males.

- gID:

  vector with IDs of SNP-genotyped individuals; rownames of genotype
  matrix in the exact order.

- DumPrefix:

  length-2 character vector to make dummy IDs; length-3 in case of
  hermaphrodites.

## Value

A character vector with IDs.
