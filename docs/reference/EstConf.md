# Confidence Probabilities

Estimate confidence probabilities ('backward') and assignment error
rates ('forward') per category (genotyped/dummy) by repeatedly
simulating genotype data from a reference pedigree using
[`SimGeno`](https://jiscah.github.io/reference/SimGeno.md),
reconstruction a pedigree from this using
[`sequoia`](https://jiscah.github.io/reference/sequoia.md), and counting
the number of mismatches using
[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md).

## Usage

``` r
EstConf(
  Pedigree = NULL,
  LifeHistData = NULL,
  args.sim = list(nSnp = 400, SnpError = 0.001, ParMis = c(0.4, 0.4)),
  args.seq = list(Module = "ped", Err = 0.001, Tassign = 0.5, CalcLLR = FALSE),
  nSim = 10,
  nCores = 1,
  quiet = TRUE
)
```

## Arguments

- Pedigree:

  reference pedigree from which to simulate, dataframe with columns
  id-dam-sire. Additional columns are ignored.

- LifeHistData:

  dataframe with id, sex (1=female, 2=male, 3=unknown), birth year, and
  optionally BY.min - BY.max - YearLast.

- args.sim:

  list of arguments to pass to
  [`SimGeno`](https://jiscah.github.io/reference/SimGeno.md), such as
  `nSnp` (number of SNPs), `SnpError` (genotyping error rate) and
  `ParMis` (proportion of non-genotyped parents). Set to `NULL` to use
  all default values.

- args.seq:

  list of arguments to pass to
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md), such as
  `Module` ('par' or 'ped'), `Err` (assumed genotyping error rate), and
  `Complex`. May include (part of) `SeqList`, a list of sequoia output
  (i.e. as a list-within-a-list). Set to `NULL` to use all default
  values.

- nSim:

  number of rounds of simulate - reconstruct - compare to perform, i.e.
  number of simulated datasets.

- nCores:

  number of computer cores to use. If `>1`, package parallel is used.
  Set to NULL to use all but one of the available cores, as detected by
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  (using all cores tends to freeze up your computer). With large
  datasets, the amount of computer memory may be the limiting factor for
  the number of cores you can use.

- quiet:

  suppress messages. `TRUE` runs `SimGeno` and `sequoia` quietly,
  `'very'` also suppresses other messages and the simulation counter
  when `nCores=1` (there is no simulation counter when `nCores>1`).

## Value

A list, with elements:

- ConfProb:

  See below

- PedErrors:

  See below

- Pedigree.reference:

  the pedigree from which data was simulated

- LifeHistData:

- Pedigree.inferred:

  a list with for each simulation the inferred pedigree based on the
  simulated data

- SimSNPd:

  a list with for each simulation the IDs of the individuals simulated
  to have been genotyped

- PedComp.fwd:

  array with `Counts` from the 'forward' `PedCompare`, from which
  `PedErrors` is calculated

- RunParams:

  a list with the call to `EstConf` as a semi-nested list (args.sim,
  args.seq, nSim, nCores), as well as the default parameter values for
  `SimGeno` and `sequoia`.

- RunTime:

  `sequoia` runtime per simulation in seconds, as measured by
  [`system.time`](https://rdrr.io/r/base/system.time.html)`()['elapsed']`.

Dataframe `ConfProb` has 7 columns:

- id.cat, dam.cat, sire.cat:

  Category of the focal individual, dam, and sire, in the pedigree
  inferred based on the simulated data. Coded as G=genotyped, D=dummy,
  X=none

- dam.conf:

  Probability that the dam is correct, given the categories of the
  assigned dam and sire (ignoring whether or not the sire is correct)

- sire.conf:

  as `dam.conf`, for the sire

- pair.conf:

  Probability that both dam and sire are correct, given their categories

- N:

  Number of individuals per category-combination, across all `nSim`
  simulations

Array `PedErrors` has three dimensions:

- class:

  - `FalseNeg`(atives): could have been assigned but was not
    (individual + parent both genotyped or dummifiable; P1only in
    `PedCompare`).

  - `FalsePos`(itives): no parent in reference pedigree, but one was
    assigned based on the simulated data (P2only)

  - `Mismatch`: different parents between the pedigrees

- cat:

  Category of individual + parent, as a two-letter code where the first
  letter indicates the focal individual and the second the parent;
  G=Genotyped, D=Dummy, T=Total

- parent:

  dam or sire

## Details

The confidence probability is taken as the number of correct (matching)
assignments, divided by all assignments made in the *observed*
(inferred-from-simulated) pedigree. In contrast, the false negative &
false positive assignment rates are proportions of the number of parents
in the *true* (reference) pedigree. Each rate is calculated separately
for dams & sires, and separately for each category
(**G**enotyped/**D**ummy(fiable)/**X** (none)) of individual, parent and
co-parent.

This function does not know which individuals in the actual `Pedigree`
are genotyped, so the confidence probabilities need to be added to the
`Pedigree` as shown in the example at the bottom.

A confidence of \\1\\ means all assignments on simulated data were
correct for that category-combination. It should be interpreted as (and
perhaps modified to) \\\> 1 - 1/N\\, where sample size `N` is given in
the last column of the `ConfProb` and `PedErrors` dataframes in the
output. The same applies for a false negative/positive rate of \\0\\
(i.e. to be interpreted as \\\< 1/N\\).

## Assumptions

Because the actual true pedigree is (typically) unknown, the provided
reference pedigree is used as a stand-in and assumed to be the true
pedigree, with unrelated founders. It is also assumed that the
probability to be genotyped is equal for all parents; in each round, a
new random set of parents (proportion set by `ParMis`) is mimicked to be
non-genotyped. In addition, SNPs are assumed to segregate independently.

An experimental version offering more fine-grained control is available
at https://github.com/JiscaH/sequoiaExtra .

## Object size

The size in Kb of the returned list can become pretty big, as each of
the inferred pedigrees is included. When running `EstConf` many times
for a range of parameter values, it may be prudent to save the required
summary statistics for each run rather than the full output.

## Errors

If you have a large pedigree and try to run this function on multiple
cores, you may run into "Cannot allocate vector of size ..." errors or
even unexpected crashes: there is not enough computer memory for each
separate run. Try reducing \`nCores\`.

## See also

[`SimGeno`](https://jiscah.github.io/reference/SimGeno.md)`, `[`sequoia`](https://jiscah.github.io/reference/sequoia.md)`, `[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md).

## Examples

``` r
# estimate proportion of parents that are genotyped (= 1 - ParMis)
prop_parents_genotyped <- c(
  dam = mean(unique(SeqOUT_griffin$Pedigree$dam) %in% rownames(Geno_griffin)),
sire = mean(unique(SeqOUT_griffin$Pedigree$sire) %in% rownames(Geno_griffin))
)

# Example for parentage assignment only
conf_grif <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
               LifeHistData = SeqOUT_griffin$LifeHist,
               args.sim = list(nSnp = 150,   # no. in actual data, or what-if
                               SnpError = 5e-3,  # best estimate, or what-if
                               CallRate=0.9,     # from SnpStats()
                               ParMis=c(0.28, 0.22)),  # calc'd above
               args.seq = list(Err=5e-3, Module="par"),  # as in real run
               nSim = 1,   # try-out, proper run >=20 (10 if huge pedigree)
               nCores=1)
#> ℹ Simulating parentage assignment only ...
#> i= 1      15:56:17 

# parent-pair confidence, per category (Genotyped/Dummy/None)
conf_grif$ConfProb
#>   id.cat dam.cat sire.cat dam.conf sire.conf pair.conf  N
#> 1      G       G        G        1         1         1 65
#> 2      G       G        X        1        NA        NA 15
#> 3      G       X        G       NA         1        NA 12
#> 4      G       X        X       NA        NA        NA 48

# Proportion of true parents that was correctly assigned
1 - apply(conf_grif$PedErrors, MARGIN=c('cat','parent'), FUN=sum, na.rm=TRUE)
#>     parent
#> cat        dam      sire
#>   GG 0.9411765 0.9390244
#>   GD 0.0000000 0.0000000
#>   GT 0.7619048 0.7857143
#>   DG 0.0000000 0.0000000
#>   DD 0.0000000 0.0000000
#>   DT 0.0000000 0.0000000
#>   TT 0.6557377 0.6637931

# add columns with confidence probabilities to pedigree
# first add columns with category (G/D/X)
Ped.withConf <- getAssignCat(Pedigree = SeqOUT_griffin$Pedigree,
                             SNPd = SeqOUT_griffin$PedigreePar$id)
Ped.withConf <- merge(Ped.withConf, conf_grif$ConfProb, all.x=TRUE,
                      sort=FALSE)  # (note: merge() messes up column order)
head(Ped.withConf[Ped.withConf$dam.cat=="G", ])
#>    id.cat dam.cat sire.cat          id         dam        sire LLRdam LLRsire
#> 30      G       G        G i110_2006_M i061_2004_F i073_2004_M   5.67    6.34
#> 31      G       G        G i025_2002_M i014_2001_F i018_2001_M   8.78    8.22
#> 32      G       G        G i136_2007_M i095_2005_F i089_2005_M   5.89    4.77
#> 33      G       G        G i145_2008_F i132_2007_F i127_2007_M   9.23    8.02
#> 34      G       G        G i079_2004_F i041_2003_F i046_2003_M  12.66    7.90
#> 35      G       G        G i080_2004_M i041_2003_F i039_2002_M  11.47    5.72
#>    LLRpair OHdam OHsire MEpair dam.conf sire.conf pair.conf  N
#> 30   15.51     0      0      0        1         1         1 65
#> 31   16.13     0      0      0        1         1         1 65
#> 32   14.83     0      0      0        1         1         1 65
#> 33   16.29     0      0      0        1         1         1 65
#> 34   17.56     0      0      0        1         1         1 65
#> 35   17.07     0      0      0        1         1         1 65

# save output summary
if (FALSE) { # \dontrun{
conf_griff[['Note']] <- 'You could add a note'
saveRDS(conf_grif[c('ConfProb','PedComp.fwd','RunParams','RunTime','Note')],
   file = 'conf_200SNPs_Err005_Callrate80.RDS')
} # }

## overall assignment rate (AR), error rate (ER) & runtime
AR_max <- sum(!is.na(Ped_griffin$dam)) + sum(!is.na(Ped_griffin$sire))
ER_max <- 2*nrow(Ped_griffin)
PCT <- conf_grif$PedComp.fwd[,'TT',,]   # Total-Total counts
list(AR = mean(apply(PCT[,'Match',],1,sum)/AR_max),  # sum over dam+sire
     ER = mean(apply(PCT[,c('Mismatch','P2only'),],1,sum)/ER_max),
     Time = mean(conf_grif$RunTime)/60)   # runtime in seconds --> minutes
#> $AR
#> [1] 0.2378788
#> 
#> $ER
#> [1] 0
#> 
#> $Time
#> [1] 0.01933333
#> 


## P(actual FS | inferred as FS) etc.
if (FALSE) { # \dontrun{
PairL <- list()
for (i in 1:length(conf_grif$Pedigree.inferred)) {  # nSim
  cat(i, "\t")
  PairL[[i]] <- ComparePairs(conf_grif$Pedigree.reference,
                             conf_grif$Pedigree.inferred[[i]],
                             GenBack=1, patmat=TRUE, ExcludeDummies = TRUE,
                             Return="Counts")
}
# P(actual relationship (Ped1) | inferred relationship (Ped2))
PairRel.prop.A <- plyr::laply(PairL, function(M)
                     sweep(M, MARGIN='Ped2', STATS=colSums(M), FUN="/"))
PairRel.prop <- apply(PairRel.prop.A, 2:3, mean, na.rm=TRUE) #avg across sims
round(PairRel.prop, 3)
# or: P(inferred relationship | actual relationship)
PairRel.prop2 <- plyr::laply(PairL, function(M)
   sweep(M, MARGIN='Ped1', STATS=rowSums(M), FUN="/"))
} # }

if (FALSE) { # \dontrun{
# confidence probability vs. sibship size
source('https://raw.githubusercontent.com/JiscaH/sequoiaExtra/main/conf_vs_sibsize.R')
conf_grif_nOff <- Conf_by_nOff(conf_grif)
conf_grif_nOff['conf',,'GD',]
conf_grif_nOff['N',,'GD',]
} # }
```
