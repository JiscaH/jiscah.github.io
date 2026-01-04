# Package index

## Input

Read, convert and check input data

- [`GenoConvert()`](https://jiscah.github.io/reference/GenoConvert.md) :
  Convert Genotype Data
- [`LHConvert()`](https://jiscah.github.io/reference/LHConvert.md) :
  Extract Sex and Birth Year from PLINK File
- [`CheckGeno()`](https://jiscah.github.io/reference/CheckGeno.md) :
  Check Genotype Matrix
- [`SnpStats()`](https://jiscah.github.io/reference/SnpStats.md) : SNP
  Summary Statistics
- [`CalcMaxMismatch()`](https://jiscah.github.io/reference/CalcMaxMismatch.md)
  : Maximum Number of Mismatches

## Simulate

Simulate genotype data assuming independent SNPs

- [`SimGeno()`](https://jiscah.github.io/reference/SimGeno.md) :
  Simulate Genotypes
- [`MkGenoErrors()`](https://jiscah.github.io/reference/MkGenoErrors.md)
  : Simulate Genotyping Errors
- [`Inherit_patterns`](https://jiscah.github.io/reference/Inherit_patterns.md)
  : Inheritance patterns

## AgePrior

See [ageprior
vignette](https://jiscah.github.io/articles/vignette_age/book/index.md)
for details

- [`MakeAgePrior()`](https://jiscah.github.io/reference/MakeAgePrior.md)
  : Age Priors
- [`PlotAgePrior()`](https://jiscah.github.io/reference/PlotAgePrior.md)
  : Plot Age Priors
- [`CalcBYprobs()`](https://jiscah.github.io/reference/CalcBYprobs.md) :
  Birth year probabilities

## Pedigree reconstruction

The core of the package, see [main
vignette](https://jiscah.github.io/articles/vignette_main/book/key-points.md)
for details

- [`sequoia()`](https://jiscah.github.io/reference/sequoia.md) :
  Pedigree Reconstruction
- [`GetMaybeRel()`](https://jiscah.github.io/reference/GetMaybeRel.md) :
  Find Putative Relatives

## Likelihood (ratios) & probabilities

- [`CalcOHLLR()`](https://jiscah.github.io/reference/CalcOHLLR.md) :
  Calculate OH and LLR for a pedigree
- [`CalcPairLL()`](https://jiscah.github.io/reference/CalcPairLL.md) :
  Calculate Likelihoods for Alternative Relationships
- [`PlotPairLL()`](https://jiscah.github.io/reference/PlotPairLL.md) :
  Plot Pair Log10-Likelihoods
- [`GetLLRAge()`](https://jiscah.github.io/reference/GetLLRAge.md) :
  LLR-age from Ageprior Matrix
- [`LLtoProb()`](https://jiscah.github.io/reference/LLtoProb.md) :
  transform log-likelihoods to probabilities
- [`CalcParentProbs()`](https://jiscah.github.io/reference/CalcParentProbs.md)
  : Calculate assignment probabilities

## Pedigree check

These functions can be applied to any pedigree, no genetic data needed

- [`EstConf()`](https://jiscah.github.io/reference/EstConf.md) :
  Confidence Probabilities
- [`GetRelM()`](https://jiscah.github.io/reference/GetRelM.md) : Matrix
  with Pairwise Relationships
- [`PlotRelPairs()`](https://jiscah.github.io/reference/PlotRelPairs.md)
  : Plot Pairwise Relationships
- [`SummarySeq()`](https://jiscah.github.io/reference/SummarySeq.md) :
  Summarise Sequoia Output or Pedigree
- [`PlotPropAssigned()`](https://jiscah.github.io/reference/PlotPropAssigned.md)
  : Plot proportion of individuals that has a parent assigned
- [`PlotSeqSum()`](https://jiscah.github.io/reference/PlotSeqSum.md) :
  Plot Summary Overview of sequoia Output
- [`getAssignCat()`](https://jiscah.github.io/reference/getAssignCat.md)
  : Assignability of Reference Pedigree
- [`getGenerations()`](https://jiscah.github.io/reference/getGenerations.md)
  : Count Generations
- [`GetAncestors()`](https://jiscah.github.io/reference/GetAncestors.md)
  : Get ancestors
- [`GetDescendants()`](https://jiscah.github.io/reference/GetDescendants.md)
  : Get descendants
- [`CalcRped()`](https://jiscah.github.io/reference/CalcRped.md) :
  Calculate Pedigree Relatedness

### Compare 2 pedigrees

- [`PedCompare()`](https://jiscah.github.io/reference/PedCompare.md) :
  Compare Two Pedigrees
- [`PlotPedComp()`](https://jiscah.github.io/reference/PlotPedComp.md) :
  Visualise PedCompare Output
- [`ComparePairs()`](https://jiscah.github.io/reference/ComparePairs.md)
  : Compare Pairwise Relationships

## Miscellaneous

- [`PedPolish()`](https://jiscah.github.io/reference/PedPolish.md) : Fix
  Pedigree
- [`ErrToM()`](https://jiscah.github.io/reference/ErrToM.md) : Generate
  Genotyping Error Matrix
- [`Err_RADseq()`](https://jiscah.github.io/reference/Err_RADseq.md) :
  Convert Genotyping Error Rates from per-allele to per-locus
- [`CountOH()`](https://jiscah.github.io/reference/CountOH.md) : Count
  opposing homozygous SNPs between pairs of individuals
- [`writeSeq()`](https://jiscah.github.io/reference/writeSeq.md) : Write
  Sequoia Output to File
- [`writeColumns()`](https://jiscah.github.io/reference/writeColumns.md)
  : Write Data to a File Column-wise
- [`FindFamilies()`](https://jiscah.github.io/reference/FindFamilies.md)
  : Assign Family IDs
- [`PedStripFID()`](https://jiscah.github.io/reference/PedStripFID.md) :
  Back-transform IDs

## Example data

- [`Ped_HSg5`](https://jiscah.github.io/reference/Ped_HSg5.md) : Example
  pedigree: 'HSg5'
- [`LH_HSg5`](https://jiscah.github.io/reference/LH_HSg5.md) : Example
  life history file: 'HSg5'
- [`SimGeno_example`](https://jiscah.github.io/reference/SimGeno_example.md)
  : Example genotype file: 'HSg5'
- [`Geno_HSg5`](https://jiscah.github.io/reference/Geno_HSg5.md) :
  Example genotype file: 'HSg5'
- [`SeqOUT_HSg5`](https://jiscah.github.io/reference/SeqOUT_HSg5.md) :
  Example output from pedigree inference: 'HSg5'
- [`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md) :
  Example pedigree: griffins
- [`LH_griffin`](https://jiscah.github.io/reference/LH_griffin.md) :
  Example life history data: griffins
- [`FieldMums_griffin`](https://jiscah.github.io/reference/FieldMums_griffin.md)
  : Example field-observed mothers: griffins
- [`Geno_griffin`](https://jiscah.github.io/reference/Geno_griffin.md) :
  Example genotype file: Griffins
- [`SeqOUT_griffin`](https://jiscah.github.io/reference/SeqOUT_griffin.md)
  : Example output from pedigree inference: griffins
- [`Conf_griffin`](https://jiscah.github.io/reference/Conf_griffin.md) :
  Example output from estimating confidence probabilities: griffins
- [`MaybeRel_griffin`](https://jiscah.github.io/reference/MaybeRel_griffin.md)
  : Example output from check for relatives: griffins
