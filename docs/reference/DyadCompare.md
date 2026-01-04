# Compare Dyads (DEPRECATED)

Count the number of half and full sibling pairs correctly and
incorrectly assigned. DEPRECATED - PLEASE USE
[`ComparePairs`](https://jiscah.github.io/reference/ComparePairs.md)

## Usage

``` r
DyadCompare(Ped1 = NULL, Ped2 = NULL, na1 = c(NA, "0"))
```

## Arguments

- Ped1:

  original pedigree, dataframe with 3 columns: id-dam-sire.

- Ped2:

  second (inferred) pedigree.

- na1:

  the value for missing parents in Ped1.

## Value

A 3x3 table with the number of pairs assigned as full siblings (FS),
half siblings (HS) or unrelated (U, including otherwise related) in the
two pedigrees, with the classification in Ped1 on rows and that in Ped2
in columns.

## See also

[`ComparePairs`](https://jiscah.github.io/reference/ComparePairs.md)
which supersedes this function;
[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
DyadCompare(Ped1=Ped_HSg5, Ped2=SeqOUT_HSg5$Pedigree)
} # }
```
