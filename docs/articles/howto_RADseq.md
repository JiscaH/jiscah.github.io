# How to: RADseq data

## Background & Caveats

`Sequoia` was initially developed for SNP array data, which after
typical quality control has a low rate of missing data (\<1%) and a very
low genotyping error rate (\<0.1%). This is quite different from RADseq
and similar sequence data, which typically have a high rate of
missingness (often \>10%) and a high rate of genotyping errors (\>5%,
see e.g. ([Bresadola et al. 2020](#ref-bresadola20))).

Because of the difference in data accuracy and completeness, the optimal
strategy for pedigree reconstruction differs. With SNP array data it is
often very obvious if and how a sample pair are related, which allows
the use of the fast hill-climbing algorithm implemented in `sequoia`,
that relies on (almost) all assignments made in previous steps being
correct. In contrast, for RADseq data the relationship between samples
is often a lot less obvious. In such situations, algorithms that fairly
randomly explore the huge number of possible pedigree configurations are
more powerful. However, this is very time consuming, and (at least some
years ago) only implemented for parentage assignment
(<https://CRAN.R-project.org/package=MasterBayes> ; archived), for
sibship clustering and parentage assignment within a single cohort
(<https://www.zsl.org/about-zsl/resources/software/colony>), or for
idealised datasets (various theoretical papers; assuming no inbreeding
and/or no missing data and/or discrete generations, etc).

Substantial power can be gained by making use of information from
neighbouring SNPs, rather than treating each SNP as independent as
`sequoia`, `Colony`, and most other programs do. However, this added
complexity is likely to come with a substantial increase in computation
time.

*Do please let me know if you a good general pedigree reconstruction
program for RADseq data, so that I can recommend it to users for which
`sequoia` does not work well!*

So, while `sequoia` is thus not an ideal tool to use with RADseq data,
it is often sufficient, and may sometimes be the only suitable tool
available.

Keeping that in mind, there are various tips & tricks to get the most
out of your data, listed below.

## Data filtering

This step is crucial; simply using all genetic data as input for
`sequoia` is almost guaranteed to give you very poor results.

Filter data based on:

- **Autosomal SNPs only**. Remove any SNPs that are heterozygous in the
  heterogametic sex (females in mammals, males in birds and some
  reptiles), but always homozygous in the other sex.

- **Minor Allele Frequency (MAF)**. SNPs with very low MAF are not
  informative, as almost all individuals will be identical homozygotes.
  Moreover, with typical RADseq error rates, any heterozygotes or
  rare-allele-homozygotes are almost as likely to be incorrect as true
  carriers of the rare allele. Therefore, inclusion of these SNPs may
  obfuscate the true pedigree signal, rather than strengthen it.

- **Missingness**. If a SNP is only scored for say 20% of individuals,
  this means that only $`0.2^2=0.04`$ or 4% of sample pairs are both
  scored, so for 96% of pairs the SNP carries no information at all.

- **Linkage Disequilibrium (LD)**. By chance, pairs of individuals will
  be more related in some regions of their genomes than you’d expect
  from the pedigree, and less related in other regions. As an extreme
  example, if you have 50 SNPs in very close LD, and 30 other SNPs
  scattered throughout the genome, the signal from those 50 SNPs can
  overwhelm the signal from the rest of the genome, and bias the
  inference of how the samples are related. It is possible to filter
  SNPs based on LD even if the genetic map of the species is unknown.

- **Genotyping error rate**. If you have some estimate of the genotyping
  error rate, for example from samples genotyped twice, consider
  removing SNPs with the very worst error rates. There is no need to
  filter extremely harshly: even if all SNPs had exactly the same
  underlying rate, you would see a Poisson distribution of the error
  count per SNP with a long tail to the right.

- **Deviation from HWE** Strong deviations from genotype proportions
  under Hardy-Weinberg Equilibrium suggest null alleles or other
  genotyping issues. Again no need to filter extremely harshly, as
  variation in this metric is expected, especially in small inbred
  populations and/or due to selection.

There is no one-size-fits-all recommendation on thresholds to use, or
the relative emphasis on the different filtering criteria – it depends
on the number of SNPs you start with, the data quality, the genome size
and chromosome number, … . It can be useful to first remove the very
worst SNPs, then remove the very worst samples, and only then do the
final filtering.

Aim for at least 100-200 SNPs for parentage assignment, and more for
sibship clustering (say 400-500). Generally, the more complicated the
pedigree, in terms of inbreeding, overlapping generations, unknown birth
years, etc., the more SNPs you will need. But, genome size and
chromosome number put an upper limit on the number of semi-independent
SNPs. Consequently, `sequoia` is unlikely to work well on inbred
populations of species with a small genome.

## Genotyping error pattern

### SNP array vs RADseq

By default, `sequoia` assumes that the chance that an allele is observed
incorrectly is independent from whether the true genotype is homozygote
or heterozygote. This is a reasonable assumption for SNP array data,
where most SNPs that are strongly affected by allelic dropout are
removed during quality control. Remaining errors largely arise from the
clustering software used to translate the signal strength for each of
the two alleles into three genotype classes. For each genotype, the
probability to be correctly observed is $`1-E`$, with for the
heterozygote an equal chance to be either homozygote (each $`E/2`$), and
for the homozygotes a much smaller chance to be observed as the other
homozygote ($`(E/2)^2`$) than as heterozygote ($`E - (E/2)^2`$).

In contrast, for RADseq the chance that a true heterozygote is observed
as homozygote is considerably higher than the reverse, due to much
higher allelic dropout rates than for SNP arrays. In ([Bresadola et al.
2020](#ref-bresadola20)), the per-allele genotyping error rate at
homozygous sites $`\epsilon_0`$ is mentioned as being typically \<1%,
while the rate at heterozygous sites $`\epsilon_1`$ is around 5-10%, and
can be much higher.

### Implementation

All relevant functions in the package allow full flexibility with
respect to the genotyping error pattern by letting you specify a vector
with three genotyping error rates:

- **hom\|hom**: the probability that a true homozygote is observed as
  the other homozygote;
- **het\|hom**: the probability that a true homozygote is observed as
  heterozygote;
- **hom\|het**: the probability that a true heterozygote is observed as
  the minor homozygote, which is equal to the probability it is observed
  as the major homozygote.

So, the only assumption made is that the two homozygotes are
interchangeable with respect to genotyping error chances.

There are two functions to help with creating these vectors,
[`ErrToM()`](https://jiscah.github.io/reference/ErrToM.md) (introduced
in version 2.0 in 2020) and
[`Err_RADseq()`](https://jiscah.github.io/reference/Err_RADseq.md)
(introduced in version 3.0 in 2025). Both can return both a length-3
vector and a 3x3 matrix; both of which can be used as input to other
functions. See the
[`?ErrToM`](https://jiscah.github.io/reference/ErrToM.md) helpfile for
further details.

## Code Examples

With a per-allele genotyping error rate at homozygous sites of
$`E_0=0.005`$ and at heterozygous sites of $`E_1=0.10`$, the error
matrix would be

``` r
Err_RADseq(E0=0.005, E1=0.10, Return='matrix')
```

    ##    obs
    ## act        0       1        2
    ##   0 0.990025 0.00995 0.000025
    ##   1 0.090000 0.82000 0.090000
    ##   2 0.000025 0.00995 0.990025

where 0, 1 and 2 are the number of copies of the reference (minor)
allele.

For comparison, the default SNP-array based pattern with the same
$`P(hom|het)=0.09`$ would be

``` r
ErrToM(0.09*2, Return='matrix')
```

    ##    obs
    ## act      0      1      2
    ##   0 0.8200 0.1719 0.0081
    ##   1 0.0900 0.8200 0.0900
    ##   2 0.0081 0.1719 0.8200

and thus has a MUCH higher $`P(het|hom)`$. Basically, with the RADseq
pattern, any observed heterozygote is probably correct
($`.82/(.82+2*.00995)= 0.976`$), but with the default pattern
$`P(het|hom)`$ is always about twice as large as $`P(hom|het)`$, and any
heterozygote is quite likely to be wrong ($`.82/(.82+2*.1719)= 0.705`$).

You can try running
[`sequoia()`](https://jiscah.github.io/reference/sequoia.md) (or any
other function) with several different error vectors, to explore the
effect it has on assignments. Hopefully, the majority of assignments is
not affected by the assumed genotyping error rate or any other input
parameter, forming a solid core of your pedigree.

``` r
Err_RAD_low <- Err_RADseq(E0=0.002, E1=0.05)
Err_RAD_high <- Err_RADseq(E0=0.01, E1=0.15)
Err_chip <- ErrToM(0.15*(1-0.15)*2, Return='vector')
```

``` r
SeqOUT_lowErr <- sequoia(GenoM, LHdata, Err=Err_RAD_low)
SeqOUT_highErr <- sequoia(GenoM, LHdata, Err=Err_RAD_high)
SeqOUT_chipErr <- sequoia(GenoM, LHdata, Err=Err_chip)
```

You can also use
[`EstConf()`](https://jiscah.github.io/reference/EstConf.md) to explore
potential consequences of the actual (i.e., simulated) genotyping error
rate being much higher/lower than assumed (during pedigree
reconstruction), or of assuming the wrong pattern:

``` r
EC_A <- EstConf(Ped_griffin, LH_griffin, 
                args.sim=list(SnpError=Err_RAD_high, nSnp=200, CallRate=0.8),
                args.seq=list(Err=Err_RAD_low),
                nSim = 10, nCores = 5)

EC_B <- EstConf(Ped_griffin, LH_griffin, 
                args.sim=list(SnpError=Err_RAD_high, nSnp=200, CallRate=0.8),
                args.seq=list(Err=Err_chip),
                nSim = 10, nCores = 5)

EC_C <- EstConf(Ped_griffin, LH_griffin, 
                args.sim=list(SnpError=Err_RAD_high, nSnp=200, CallRate=0.8),
                args.seq=list(Err=Err_RAD_high),
                nSim = 10, nCores = 5)
```

Calculating the number of correct and incorrect assignments, averaged
across replicates, and summed over dams and sires:

``` r
count_ARER <- function(EC) {
  list(AR = mean(apply(EC$PedComp.fwd[,'TT','Match',],1,sum)),
       ER = mean(apply(EC$PedComp.fwd[,'TT',c('Mismatch','P2only'),], 1, sum)))
}

count_ARER(EC_A)  # RAD pattern, presumed rate lower than simulated
# $AR
# [1] 71.3
# $ER
# [1] 27.2

count_ARER(EC_B)  # SNP chip pattern, similar/higher rate presumed than sim
# $AR
# [1] 250.8
# $ER
# [1] 47.8

count_ARER(EC_C)  # presumed rate and pattern identical to simulated
# $AR
# [1] 266.3
# $ER
# [1] 8.2
```

Of course, in reality you likely will not know the actual genotyping
error rates very accurately (if at all), but trying out different values
in [`sequoia()`](https://jiscah.github.io/reference/sequoia.md),
[`EstConf()`](https://jiscah.github.io/reference/EstConf.md),
[`GetMaybeRel()`](https://jiscah.github.io/reference/GetMaybeRel.md),
[`CalcOHLLR()`](https://jiscah.github.io/reference/CalcOHLLR.md) and/or
[`CalcPairLL()`](https://jiscah.github.io/reference/CalcPairLL.md) will
get you a feeling of what might and might not be possible with your
dataset.

## Literature

Bresadola, L, V Link, C A Buerkle, C Lexer, and D Wegmann. 2020.
“Estimating and Accounting for Genotyping Errors in RAD-Seq
Experiments.” *Molecular Ecology Resources* 20 (4): 856–70.
