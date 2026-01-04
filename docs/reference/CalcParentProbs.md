# Calculate assignment probabilities

For each assigned offspring-parent pair, calculate the probability they
are parent-offspring vs otherwise related. Probabilities are scaled to
sum to one across all possible\* relationships between the pair or trio;
see Details.

## Usage

``` r
CalcParentProbs(Pedigree = NULL, GenoM = NULL, quiet = FALSE, nCores = 1, ...)
```

## Arguments

- Pedigree:

  dataframe with columns id-dam-sire. By default, any non-genotyped
  individuals are 'dummified'; use `Module='par'` to ignore them.

- GenoM:

  numeric matrix with genotype data: One row per individual, one column
  per SNP, coded as 0, 1, 2, missing values as a negative number or NA.
  You can reformat data with
  [`GenoConvert`](https://jiscah.github.io/reference/GenoConvert.md), or
  use other packages to get it into a genlight object and then use
  `as.matrix`.

- quiet:

  logical, suppress messages. No progress is printed when \>1 core is
  used.

- nCores:

  number of computer cores to use. If `2` or `4`, package parallel is
  used (other values are not applicable).

- ...:

  Additional arguments passed to
  [`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md), such
  as the genotyping error rate `Err`, age information in `LifeHistData`
  and `AgePrior`, or `InclDup` to include the probability that the two
  samples are duplicates.

## Value

the `Pedigree` dataframe with the three applicable columns renamed to
id-dam-sire, and 7 additional columns:

- Probdam:

  Probability that individual in dam column is the maternal parent,
  rather than otherwise related (LL(PO)/sum(LL))

- Probsire:

  Analogous for sire

- Probpair:

  Probability for id-dam-sire trio. Approximated as the minimum of dam
  conditional on sire and sire conditional on dam, thus not including
  e.g. both being siblings (those other configurations are considered by
  sequoia during pedigree reconstruction, but can (currently) not be
  accessed directly)

- dam_alt, sire_alt:

  Most likely alternative (not PO) relationship between id-dam and
  id-sire, respectively

- Probdam_alt, Probsire_alt:

  Probability of most likely alternative relationship

## Details

The returned probabilities are calculated from the likelihoods used
throughout the rest of this package, by scaling them to sum to one
across all possible relationships. For `Complex='simp'` these are
PO=parent-offspring, FS=full siblings, HS=half siblings,
GP=grand-parental, FA=full avuncular, HA=third degree relatives (incl
half avuncular), and U=unrelated. For `Complex='full'` there are
numerous double relationship considered (PO & HS, HS & HA, etc), making
both numerator and denominator in the scaling step less unambiguous, and
the returned probabilities an approximation.

The likelihoods are calculated by calling
[`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md) once or
twice for each id-dam and id-sire pair: once not conditioning on the
co-parent, and once conditional on the co-parent, if any. For genotyped
individuals this is done with `focal='PO'`, and for dummy individuals
with `focal='GP'`.

For relationships between a genotyped and a dummy individual, it may
only be possible to determine that the genotyped individual is a second
degree relative (GP, HS, or FA) to the dummy's offspring. This then
results in a probability of at most 0.33, even when the two are indeed
parent and offspring.

See [`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md) and
the vignettes for further details.

Note that for large pedigrees this function can be fairly slow,
especially when using
[`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md)'s
default `Module='ped'` and `Complex='full'`.

Subsetting the genotype data may give different results, as the
likelihoods and thus the probabilities depend on the allele frequencies
in the sample.

## Warning

The probabilities will be less reliable with close inbreeding and double
relationships. This function has not been tested yet with
hermaphrodites, and is unlikely to give reliable results without further
code updates.

## See also

[`CalcPairLL`](https://jiscah.github.io/reference/CalcPairLL.md),
[`LLtoProb`](https://jiscah.github.io/reference/LLtoProb.md)

## Examples

``` r
test_ped <- Ped_griffin[21:25,]
# add an incorrect sire to illustrate
test_ped$sire <- as.character(test_ped$sire)
test_ped$sire[5] <- 'i057_2003_M'
Ped_with_probs <- CalcParentProbs(test_ped, Geno_griffin)
#> ! `GenoM`/`gID` and `Pedigree` share few common individuals
#>  
#> ℹ Calculating LL for  dam_solo
#> ℹ Conditioning on pedigree with 145 individuals, 5 dams and 4 sires
#> Transferring input pedigree ...
#>  
#> ℹ Calculating LL for  sire_solo
#> ℹ Conditioning on pedigree with 145 individuals, 5 dams and 4 sires
#> Transferring input pedigree ...
#>  
#> ℹ Calculating LL for  dam_paired
#> ℹ Conditioning on pedigree with 145 individuals, 5 dams and 4 sires
#> Transferring input pedigree ...
#>  
#> ℹ Calculating LL for  sire_paired
#> ℹ Conditioning on pedigree with 145 individuals, 5 dams and 4 sires
#> Transferring input pedigree ...
print(Ped_with_probs, digits=2)
#>            id         dam        sire birthyear Probdam dam_alt Probdam_alt
#> 1 i021_2002_M i001_2001_F i020_2001_M      2002       1      FS     2.4e-06
#> 2 i022_2002_F i015_2001_F i006_2001_M      2002       1      GP     2.8e-09
#> 3 i023_2002_F i007_2001_F i019_2001_M      2002       1      FS     7.9e-07
#> 4 i024_2002_M i013_2001_F        <NA>      2002      NA    <NA>          NA
#> 5 i025_2002_M i014_2001_F i057_2003_M      2002       1      GP     2.3e-08
#>   Probsire sire_alt Probsire_alt Probpair
#> 1    1e+00       HS      2.4e-08  1.0e+00
#> 2    1e+00       GP      9.3e-08  1.0e+00
#> 3       NA     <NA>           NA       NA
#> 4       NA     <NA>           NA       NA
#> 5    6e-34       HS      2.5e-01  6.7e-61
# Any non-genotyped non-'dummifiable' individuals are automatically skipped

# To get likelihoods for 'all' relationships, not just probabilities for
# PO & (next-)most-likely:
LL_sire_single <- CalcPairLL(
  Pairs = data.frame(id1=test_ped$id,
                     id2=test_ped$sire,
                     dropPar1='both', # drop both -> id2 as single parent
                     focal='PO'),
  Pedigree = Ped_griffin,   # pedigree to condition on
  GenoM = Geno_griffin, Plot=FALSE)
#> ℹ Conditioning on pedigree with 200 individuals, 167 dams and 163 sires
#> ! Assuming columns 1 and 2 of Pairs are 'ID1' and 'ID2'
#> Transferring input pedigree ...
```
