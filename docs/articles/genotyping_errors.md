# Genotyping errors

No genotyping method is 100% error free. For some methods, such as SNP
arrays, the genotyping error rate $`E`$ is very low (say 0.01%). In
those cases, the assumed error pattern will have no noticeable effect.
For other methods, such as RADseq, the error rate is typically much
higher (1%, or even 5%), and then the assumed pattern *can* have a
noticeable effect.

## Per-locus vs per-allele error rate

The error rate $`E`$ (`Err`) throughout the `sequoia` documentation
represents the per-locus error rate; the per-allele error rate is
$`E/2`$. There appears to be no consistency in the literature on whether
the per-locus or per-allele error rate is used, so it is worth making
sure it is defined consistently throughout your analysis pipeline.

## Error patterns

With ‘error pattern’ I mean the relative chances of a true homozygote
**AA** being erroneously scored as a heterozygote **Aa** versus the
other homozygote **aa**, versus a true **Aa** being observed as **AA**
or **aa**.

The pattern of genotyping errors will depend on the genotyping method
and the software used to call SNPs. For SNP arrays, the clustering based
on the intensity of the **A** vs **a** probe results in the chance of a
true **AA** being scored as **Aa** being about as likely as the reverse.
In contrast, for sequencing-based methods scoring a true **Aa**
erroneously as a homozygote is much more likely than the reverse,
especially with low coverage.

Because there is thus no one-size-fits-all error model, `sequoia` offers
full flexibility in specifying the error pattern. This can be done
either

- via a length 3 vector, if the two alleles are interchangeable with
  respect to error rate (so that e.g. \$P(**aa**\|**Aa**) ==
  P(**AA**\|**Aa**)), or
- via a 3x3 matrix, which only assumes that the heterozygotes **Aa** and
  **aA** are interchangeable.

## Default model

The default error model is based on SNP array genotyping, and from
version 2.9 onwards is defined as:

|     | aa            | Aa           | AA            |
|:----|:--------------|:-------------|:--------------|
| aa  | $`(1-E/2)^2`$ | $`E(1-E/2)`$ | $`(E/2)^2`$   |
| Aa  | $`E/2`$       | $`1-E`$      | $`E/2`$       |
| AA  | $`(E/2)^2`$   | $`E(1-E/2)`$ | $`(1-E/2)^2`$ |

Probability of observed genotype (columns) conditional on actual
genotype (rows) and per-locus error rate E

Which for an error rate of 5% translates to:

``` r
sequoia::ErrToM(Err=0.05) %>% round(4)
```

    ##    obs
    ## act      0      1      2
    ##   0 0.9500 0.0494 0.0006
    ##   1 0.0250 0.9500 0.0250
    ##   2 0.0006 0.0494 0.9500

where 0, 1 or 2 refers to the number of copies of the reference allele
**A**.

Under this model, the probability that the genotype is observed
correctly is $`1-E`$, irrespective of the true genotype. The probability
that a genotype is correct, *does* differ between observed heterozygotes
and observed homozygotes:

``` r
ErM <- sequoia::ErrToM(Err=0.05)
round( 0.95/colSums(ErM), 3)
```

    ##     0     1     2 
    ## 0.974 0.906 0.974

So an observed heterozygote has a 90.6% chance to be correct, while a
observed homozygote a 97.4% chance.

### Error model for sequence-based methods

For RADseq and similar methods, any observed homozygote may be a
heterozygote with relatively high probability, especially with low
coverage, while any heterozygote is much less likely to be incorrect.

We can specify this for example as follows \[^3\] (numbers are
completely random!!):

``` r
ErrV_RAD <- c(0.002^2, # hom|hom: true hom observed as other hom
              0.002,   # het|hom: true hom observed as het
              0.05)    # hom|het: true het observed as hom 
```

for an allelic dropout rate of 5%, and an error rate of 0.2% per allele,
giving a hom-to-other-hom error rate of $`0.001^2 = 1e-6`$. This can be
turned into a error matrix for double checking:

``` r
ErM_RAD <- sequoia::ErrToM(Err = ErrV_RAD, Return='matrix')
ErM_RAD
```

    ##    obs
    ## act        0     1        2
    ##   0 0.997996 0.002 0.000004
    ##   1 0.050000 0.900 0.050000
    ##   2 0.000004 0.002 0.997996

``` r
round( diag(ErM_RAD) / colSums(ErM_RAD), 3)
```

    ##     0     1     2 
    ## 0.952 0.996 0.952

so now we have an error model where any observed heterozygote has a 0.4%
probability to be incorrect, and any observed homozygote a 4.8%
probability.

## Specify error matrix

You can pass a 3x3 matrix to `Err`, with the only requirement that rows
must sum to 1 (each actual genotype is observed as either **aa**,
**Aa**, or **AA**). You can use this if for example the homozygote for
the reference allele may be more likely to be mis-typed than the
homozygote for the alternative allele, the heterozygote may be mis-typed
with different probabilities for the two homozygotes, or whatever else
is appropriate for your dataset.

For example, you may have the following counts from genotyping some
samples several times:

``` r
ErrCount <- matrix(c(14, 3, 5,
                     4, 28, 9,
                     0, 4, 72),
                   nrow=3, byrow=TRUE)
ErrCount
```

    ##      [,1] [,2] [,3]
    ## [1,]   14    3    5
    ## [2,]    4   28    9
    ## [3,]    0    4   72

A problem here is the cell with $`0`$: by leaving this as is, it would
imply that true **AA** can never be observed as **aa**, even if we would
genotype thousands of samples. That seems unlikely, so as a quick fix
let’s pretend that if we had genotyped twice as many, we would have
found 1 error:

``` r
ErrCount[3,1] <- 0.5
```

Then divide by rowsums, e.g. using
[`sweep()`](https://rdrr.io/r/base/sweep.html):

``` r
ErrProbs <- sweep(ErrCount, MARGIN=1, STATS=rowSums(ErrCount), FUN="/")
ErrProbs %>% round(3)
```

    ##       [,1]  [,2]  [,3]
    ## [1,] 0.636 0.136 0.227
    ## [2,] 0.098 0.683 0.220
    ## [3,] 0.007 0.052 0.941

and then you can use this as input for the
[`sequoia()`](https://jiscah.github.io/reference/sequoia.md) analysis:

``` r
SeqOUT <- sequoia(GenoData, LifeHistData,
                  Err = ErrProbs)
```

In reality you probably want something a bit more sophisticated,
e.g. relying on the known error pattern for your genotyping platform.

### Do not use 0

Using a value of $`0`$ in the error matrix is not advised, as it puts an
infinite amount of weight on the locus if it does happen after all (see
also [Cromwell’s rule](https://en.wikipedia.org/wiki/Cromwell%27s_rule).

## Specify a function (`ErrFlavour`)

If you want to run
[`sequoia()`](https://jiscah.github.io/reference/sequoia.md) with
several presumed error rates but with the same pattern, you can create a
function that takes the value(s) as input, and returns either a length 3
vector (hom\|hom, het\|hom, hom\|het), or a 3x3 matrix with rows summing
to 1. If this function takes a single value as input, you can pass it
directly to [`sequoia()`](https://jiscah.github.io/reference/sequoia.md)
or various other functions. Those will internally call
[`ErrToM()`](https://jiscah.github.io/reference/ErrToM.md), which takes
care of the conversion to a 3x3 matrix.

For example,

``` r
RAD_flavour <- function(E)  c(E^2, E, .05+E)
test <- sequoia(Geno_griffin, LH_griffin, Module='dup', quiet=TRUE,
                Err=0.01, ErrFlavour=RAD_flavour)
```

you can inspect the error matrix with

``` r
test$ErrM
```

    ##    obs
    ## act      0    1      2
    ##   0 0.9899 0.01 0.0001
    ##   1 0.0600 0.88 0.0600
    ##   2 0.0001 0.01 0.9899

### Per-allele error rate

If you would want to keep the default error pattern, but specify `Err`
as the per-allele rather than per-locus error rate, you may do this as
follows:

``` r
ErFunc_per_allele <- function(Ea) {
  c(Ea^2, 2*Ea - Ea^2, Ea)
}

ErFunc_per_allele(0.01)
```

    ## [1] 0.0001 0.0199 0.0100

``` r
ErrToM(Err=0.01, flavour = ErFunc_per_allele, Return = 'matrix')
```

    ##    obs
    ## act      0      1      2
    ##   0 0.9800 0.0199 0.0001
    ##   1 0.0100 0.9800 0.0100
    ##   2 0.0001 0.0199 0.9800

``` r
# and for comparison:
ErrToM(Err=0.02, Return = 'vector')
```

    ## [1] 0.0001 0.0199 0.0100

### Inbuilt functions

Currently the only inbuilt functions are ones that return the error
pattern in earlier versions of `sequoia`. These are described in
[`?ErrToM`](https://jiscah.github.io/reference/ErrToM.md), and can be
used by specifying e.g. `ErrFlavour = 'version1.3'`.

If you have a good reference for error patterns for RADseq etc, please
let me know as I’d be more than happy to build those in.

## Variation between SNPs

Currently a constant error rate across all SNPs is assumed. It is
unclear whether setting this presumed rate to the arithmetic mean rate
(the regular average) across SNPs gives the best performance for real
world data, or if e.g. the harmonic mean or median may be more
appropriate.

Anecdotal evidence suggests that using a somewhat higher presumed error
rate increases performance (more correct assignments and/or fewer
incorrect assignments), but this has not been thoroughly explored. It
may partly be due to confusion between the per-locus and per-allele
error rate, but also due to the very skewed nature of most distributions
(many alleles with low error rate, and some with high error rate),
resulting in a arithmetic mean that is lower than the median.

## Effect of `Err` on pedigree reconstruction

For each individual at each SNP the probability that the actual genotype
is 0, 1 or 2 copies of the reference allele is estimated, based on among
others its own observed genotype and the genotyping error rate $`E`$. If
$`E`$ is close to zero, the actual genotype is assumed to be known with
high certainty, and the likelihood differences between various
relationships are large. The higher $`E`$, the fuzzier the differences
between the relationships become, and the more challenging pedigree
reconstruction is. This can be explored by running `CalcPairLL` a few
times on the same pedigree with the same genetic data, only varying
`Err`.

All details of calculation of likelihoods incorporating genotyping
errors are described in Huisman ([2017](#ref-huisman17)).

Noteworthy is that the probabilistic estimate of an individual’s true
genotype depends on its own observed genotype and its parents’ estimated
true genotypes, but *not on the genotypes of its offspring*. Doing so
would potentially create an unstable seesaw situation, where assignment
of a parent is both affected by and does affect that parent’s estimated
genotype.

For dummy individuals, which do not have an observed genotype, the
offspring genotypes do necessarily affect its estimated genotype. To
prevent the algorithm from seesawing, marginal probabilities are
calculated, excluding grandparents and/or one, two or all offspring from
the calculation, as appropriate for a particular relationship likelihood
(most importantly avoiding implicitly conditioning on oneself).

The main reason those marginal probabilities are not implemented for
genotyped individuals, is because of their relatively large
computational cost, with negligible effect on the pedigree outcome for
moderately high quality SNP data (error rate \< 1% and low
missingness).[¹](#fn1)

## Try out different values

Due to both the variation in error rate between SNPs and how genotyping
errors are implemented, it is usually a good idea to run
[`sequoia()`](https://jiscah.github.io/reference/sequoia.md) with a few
different presumed genotyping error rates. This allows you to make sure
the results are robust against slight mis-specifications of this
parameter, and/or to optimise the false negative and false positive
assignment rates.

### `Err = 0`

[`sequoia()`](https://jiscah.github.io/reference/sequoia.md) is not
guaranteed to work with an assumed error rate of $`0`$, especially not
for more complex pedigrees. The problem is that a new assignment
considers only the direct vicinity of the putative parent-offspring
pair, not the full pedigree. An assignment will have cascading effects
on the probabilities of the true, underlying genotypes of dummy or
partially genotyped individuals downstream in the pedigree, which may
result in impossibilities when `Err=0` but genotyping errors are
present.

------------------------------------------------------------------------

### References

Huisman, Jisca. 2017. “Pedigree Reconstruction from SNP Data: Parentage
Assignment, Sibship Clustering and Beyond.” *Molecular Ecology
Resources* 17 (5): 1009–24.

------------------------------------------------------------------------

1.  If you are interested in a change in this, e.g. as an optional
    feature, please do contact me to discuss possibilities. I am fairly
    sure it would improve performance for low quality SNP data.
