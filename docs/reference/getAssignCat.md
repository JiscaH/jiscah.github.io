# Assignability of Reference Pedigree

Identify which individuals are SNP genotyped (G), and which can
potentially be substituted by a dummy individual ('dummifiable', D).

## Usage

``` r
getAssignCat(Pedigree, SNPd, minSibSize = "1sib1GP")
```

## Arguments

- Pedigree:

  dataframe with columns id-dam-sire. Reference pedigree.

- SNPd:

  character vector with ids of genotyped individuals.

- minSibSize:

  minimum requirements to be considered dummifiable is 1 genotyped
  offspring, and

  - '1sib1GP': at least 1 grandparent (G or D) or 1 more offspring (G or
    D); these are potentially assignable by
    [`sequoia`](https://jiscah.github.io/reference/sequoia.md)

  - '2sib': at least 1 more offspring (i.e. 2 siblings). Old default for
    [`PedCompare`](https://jiscah.github.io/reference/PedCompare.md).

  .

## Value

The `Pedigree` dataframe with 3 additional columns, `id.cat`, `dam.cat`
and `sire.cat`, with coding similar to that used by
[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md):

- G:

  Genotyped

- D:

  Dummy or 'dummifiable'

- X:

  Not genotyped and not dummifiable

## Details

Non-genotyped individuals can potentially be substituted by a dummy
during pedigree reconstruction by
[`sequoia`](https://jiscah.github.io/reference/sequoia.md) when they
have at least one genotyped offspring, and either one additional
offspring (genotyped or dummy) or an genotyped/dummy parent (i.e. a
grandparent to the genotyped offspring).

Note that this is the bare minimum requirement; e.g. grandparents are
often indistinguishable from full avuncular (see
[`sequoia`](https://jiscah.github.io/reference/sequoia.md) and vignette
for details). G-G parent-offspring pairs are only assignable if there is
age information, or information from the surrounding pedigree, to tell
which of the two is the parent.

It is assumed that all individuals in `SNPd` have been genotyped for a
sufficient number of SNPs. To identify samples with a too-low call rate,
use [`CheckGeno`](https://jiscah.github.io/reference/CheckGeno.md). To
calculate the call rate for all samples, see the examples below.

## Examples

``` r
PedA <- getAssignCat(Ped_HSg5, rownames(SimGeno_example))
tail(PedA)
#>          id    dam   sire id.cat dam.cat sire.cat
#> 995  b05187 a04045 b04098      X       X        X
#> 996  a05188 a04045 b04098      X       X        X
#> 997  a05189 a04006 b04177      X       X        X
#> 998  b05190 a04006 b04177      X       X        X
#> 999  b05191 a04006 b04177      X       X        X
#> 1000 b05192 a04006 b04177      X       X        X
table(PedA$dam.cat, PedA$sire.cat, useNA="ifany")
#>    
#>       D   G   X
#>   D   4  52   0
#>   G   8 232  24
#>   X   0  64 616

# calculate call rate
if (FALSE) { # \dontrun{
CallRates <- apply(MyGenotypes, MARGIN=1,
                   FUN = function(x) sum(x!=-9)) / ncol(MyGenotypes)
hist(CallRates, breaks=50, col="grey")
GoodSamples <- rownames(MyGenotypes)[ CallRates > 0.8]
# threshold depends on total number of SNPs, genotyping errors, proportion
# of candidate parents that are SNPd (sibship clustering is more prone to
# false positives).
PedA <- getAssignCat(MyOldPedigree, rownames(GoodSamples))
} # }
```
