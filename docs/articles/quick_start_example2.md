# Quick start example 2: Real data

------------------------------------------------------------------------

First get a subset of SNPs which are as informative (high MAF, high call
rate), reliable (low error rate), and independent (low LD) as possible.
Ideally around 400 – 700 for full pedigree reconstruction; fewer are
necessary for only parentage assignment, while more may be necessary
with high levels of polygamy or inbreeding.

One possible tool to perform such subsetting of SNPs is
[PLINK](https://www.cog-genomics.org/plink2) using something like this
(on the command line), with thresholds adapted to your particular
dataset.

``` bash
plink --file mydata --geno 0.1 --maf 0.3 --indep 50 5 2
plink --file mydata --extract plink.prune.in --recodeA --out inputfile_for_sequoia
```

In addition you need a dataframe with the sex and birth/hatching year of
as many individuals as possible, in arbitrary order.

Then, in R start with something similar to this, adapting the code to
match your file names and guestimated genotyping error rate.

``` r
# read in genotype data if already coded as 0/1/2, with missing=-9:
Geno <- as.matrix(read.csv("mydata.csv", header=FALSE, row.names=1))
CheckGeno(Geno)
# read in many other input formats (not .vcf (yet)):
Geno <- GenoConvert(InFile = "mydata.ped", InFormat="ped")
#
# read in lifehistory data: ID-Sex-birthyear, column names ignored
# optional: minimum & maximum birth year, when not exactly known
LH <- read.table("LifeHistoryData.txt", header=T)
#
# duplicate check & parentage assignment (takes few minutes)
# (maximum number of sibship-clustering iterations = 0)
ParOUT <- sequoia(GenoM = Geno,  LifeHistData = LH_HSg5,
                  Module="par", Err=0.005,
                  quiet = FALSE, Plot = TRUE)
#
# inspect duplicates (intentional or accidental)
ParOUT$DupGenotype
#
# compare assigned parents to field pedigree (check column order!)
FieldPed <- read.table("FieldPed.txt", header=T)
PC.par <- PedCompare(Ped1 = FieldPed[, c("id", "dam", "sire")],
                     Ped2 = ParOUT$PedigreePar)
PC.par$Counts["TT",,]
#
# calculate Mendelian errors per SNP (works also w field pedigree)
stats <- SnpStats(Geno, ParOUT$PedigreePar)
MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1-stats[,"AF"])
#
# ..........................................................
# polish dataset: remove one indiv. from each duplicate pair
# & drop low call rate samples
# & drop SNPs with high error rate and/or low MAF
Geno2 <- Geno[!rownames(Geno),] 
Geno2 <- Geno2[, -which(stats[,"Err.hat"]>0.05 | MAF < 0.1)]
#
Indiv.Mis <- apply(Geno2, 1, function(x) sum(x == -9)) / ncol(Geno2)
Geno2 <- Geno2[Indiv.Mis < 0.2, ]
# check histograms for sensible thresholds, iterate if necessary
#
# run full pedigree reconstruction (may take up to a few hours)
# including re-run of parentage assignment
SeqOUT <- sequoia(GenoM = Geno2,
                  LifeHistData = LH_HSg5,
                  Module = "ped",
                  Err = 0.001)
#
# inspect assigned parents, proportion dummy parents, etc.
SummarySeq(SeqOUT)
# (see Example 1 for saving results)
```
