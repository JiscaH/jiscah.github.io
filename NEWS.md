
# sequoia 2.4.1

### New features & major changes
- R markdown file to create (or serve as a first draft of) a summary report of input and output, called via `sequoia_report()`
- `PedCompare()` parameter `minSibSize` (minimum sibship size for which the non-genotyped parent is considered 'dummyfiable') changed default value from `2sib` (2 genotyped siblings) to `1sib1GP`. This reflects the increased success of reconstructing grandoffspring-grandparent pairs in newest version, which make the output of `PedCompare(,minSibSize='2sib')` confusing. This also affects `EstConf()`. 
- `GetRelM()` now allows input of Pedigree AND Pairs. This a.o. allows `PlotRelPairs()` with both inferred pedigree plus `GetMaybeRel()` output.  



### Minor changes
- switch from library xlsx to openxlsx in `writeSeq()` (easier to install)
- adds purl=FALSE to vignette chunks & include `SeqOUT_HSg5` for quicker vignette compilation
- include additional griffin example data: `Geno_griffin`, `Conf_griffin` (output from `EstConf()`) and `MaybeRel_griffin` (output from `GetMaybeRel()`), as well as the script used to create these (`mk_griffin_data.R`)
- many examples are rewritten, for clarification & to speed up package check
- the plotting function of `SummarySeq()` was internal and is now exported (`PlotSeqSum()`)
- output from `sequoia()` now includes `args.AP`. 
- in `SummarySeq()`, if `Pedigree` is provided rather than `SeqList` and `SNPd=NULL`, all individuals are categorised as `Observed` (was: `Genotyped`).
- `PedPolish()` now has arguments to specify whether to drop extra columns (besides id-dam-sire) and whether to keep rows with non-unique or NA ids. 
- `SnpStats()` now includes HWE tests


### Bug fixes
- fixes several memory issues identified by CRAN
- fixes various bugs in `SimGeno()` for non-autosomal inheritance; this option is still experimental and non-autosomal SNPs are not supported by the pedigree reconstruction.
- fixes `cannot allocate vector of size ...` issue in `GetRelM()` with very large pedigrees, which affected `MakeAgePrior()` and thereby `sequoia()`. 
- `EstConf()` `$ConfProb` used the wrong denominator, namely the number of parents in the reference pedigree rather than in the inferred pedigree. 
- fixes bug resulting in LL(FA)=777 ('impossible') for pairs with all 4 parents unknown. Likely only affected `CalcPairLL()`, as LL(HS)=LL(GP)=LL(FA) for these pairs. Origin time unknown. 
- fixes bug in `GenoConvert()` when `Informat='single'`. 


# sequoia 2.3.5
- fixes bug: OH count was always zero when there was no co-parent


# sequoia 2.3.4
- fixes bug in `CalcPairLL()` HS likelihood when conditioning on pedigree was incorrect. No/minimal effect on pedigree reconstruction.
- fixes bug in `DuplicateCheck()` (always automatically called by `sequoia()`) that on very rare occasions caused R to crash


# sequoia 2.3.3
fixes minor bugs identified by CRAN valgrind and gcc-ASAN


# sequoia 2.3.2
minor edits to vignette to comply with CRAN precheck


# sequoia 2.3.1

### New features & major changes
- Hermaphrodites: dummy individuals with offspring as both dam and as sire now 
have prefix 'H'; closer links between the two 'clonal' sibship parts during pedigree reconstruction for improved performance
- Assignment of sibship grandparents moved to before check & assignment of additional parents; this proved to increase correct assignments without increasing incorrect assignments. 
- new function `CalcBYprobs()` to estimate the probability that individual i is born in year y. 

### Minor changes
- various edits to Fortran code improving general performance. 

### Bug fixes
- fixed inconsistent rounding in `EstConf()` output
- bug in `GenoConvert()` regarding InData vs InFile


# sequoia 2.3.0

### Major bug fixes
- now possible to run with assumed genotyping error rate `Err=0`
- further improvements when a large proportion of birth years are missing
- various bugs related to hermaphrodites

### Reduced computational time
- Sibship clustering etc. now done from oldest to youngest individual instead of in order of occurrence in genotype file, giving quicker convergence to final high likelihood pedigree.



# sequoia 2.2.2

### New features & major changes
- hermaphrodites: re-implemented and greatly improved sibship clustering. Specification of hermaphrodite vs diocious is now separate input parameter `Herm` instead of specified via `Complex`. New output list element `DummyClones`
- improved performance when a large proportion of birth years are missing
- new function `CalcRped()` to calculate pedigree relatedness. Uses package `kinship2`. 

### Bug fixes
Various smaller bugs have been fixed, some affecting assignment rate



# sequoia 2.1.0 

### New features & major changes
- parameter 'maxSibIter' (-9/-1/0/>0) deprecated and replaced by 'Module' (pre/dup/par/ped). 
- `sequoia()`: Option `FindMaybeRel` has been deprecated; call `GetMaybeRel()` directly instead.
- new function `CalcPairLL()`, returns likelihoods for each of the 7 considered relationships (PO, FS, HS, GP, FA, HA, U) for each specified pair of individuals
- new function `RelPlot()` for Colony-like visualisation of pairwise relationships (automatically called by `ComparePairs()`) 


### Bug fixes
- `sequoia()`: fixed error `object 'ErrM' not found` when re-using output from a previous sequoia run. Circumvent this bug in version 2.0.7 by fooling the program to think it's output from an older version: `names(ParOUT$Specs)[match("MaxMismatchOH", names(ParOUT$Specs))] <- "foo"`. 
- `sequoia()`: fixed bug causing genotyped parents to not always be monogamous when `Complx='mono'`.
- various functions: fixed error when dummy prefixes have different number of characters (`Error in data.frame(id = c(dID[s(nd[1]), 1], dID[s(nd[2]), 2]), VtoM(TMP$dumparrf,  : arguments imply differing number of rows)` )
- `GetMaybeRel()`: fixed error `(subscript) logical subscript too long` when input pedigree contains dummies
- `GetMaybeRel()`: fixed error causing likely GP pairs not to be included output
- `PedCompare()`: fixed id.dam.cat and id.sire.cat being 'NANA' instead of XD, XG or XX when Symmetrical=TRUE
- `PlotAgePrior()`: Avoid using `grDevices::hcl.colors()` in R versions <3.6, where this function is not yet available 
- `ComparePairs()`: fixed bug when `Pairs2` but not `Ped2` is specified
- `MakeAgePrior()`: fixed bug when there are no FS pairs in the input pedigree
- `MakeAgePrior()`: `MaxAgeParent` was ignored when a pedigree with overlapping generations was supplied


### Minor changes
- most examples are now set to `\donttest` instead of `\dontrun`, so that they can be run using `example()`. Note that some can be quite time consuming, especially `EstConf()`.
- new function `PlotPedComp()` to visualise `PedCompare()` output
- `SimGeno()`: deprecated input parameters (since v 1.3.1) are dropped completely
- `getAssignCat()` no longer drops additional columns from input pedigree
- Speed increase for `CalcOHLLR(, CalcLLR = FALSE)`
- More thorough input checks, which are more consistent across different functions
- `PedCompare()` output element `DummyMatch` now also include output class of matched individual's parents & offspring
- Duplicate check in `sequoia()` now only returns pairs for which LL_duplicate - max(LL_{not duplicate}) > T_filter; when call rates are low this may be a substantially shorter list than in previous versions, where all pairs with fewer than MaxMismatchDUP differences were listed. 
- `ComparePairs()` can now be called for a single pedigree, as well as to compare two pedigrees
- If the plotting window in Rstudio is too small and `Plot=TRUE`, all functions will print a message and return results as usual, instead of throwing an error and not returning results.
- `MakeAgePrior()`: more consistent implementation; is now called by `sequoia()` with only  lifehistory data of genotyped individuals. 
- `EstConf()` now also returns the full `Counts` table from `PedCompare()`; `$RunParams` now holds the evaluated input paramters, instead of e.g. `V[i]` when called from inside a loop. 


# sequoia 2.0.7 

### New features
- The genotyping error matrix (probability of observed genotype conditional on actual) is now fully customisable in all relevant functions, see help file of new function `ErrToM`. The default has changed very slightly from version 1.3. 
- new function `CalcOHLLR()` to calculate Mendelian errors and parental log-likelihood ratios for any pedigree
- new function `getAssignable()` to flag genotyped and 'dummifiable' individuals in any pedigree
- new function `ComparePairs()` to compare pairwise relationships between 2 pedigrees; replaces now-deprecated `DyadCompare`. 
- function `PedPolish()` is now user available.

### Major changes
- Deprecated option `MaxMismatch` of function `sequoia`, now calculated internally by new function `CalcMaxMismatch` based on number of SNPs, presumed genotyping error rate, and minor allele frequencies
- function `EstConf()` now also estimates confidence for parent-pairs; output has changed considerably. 
- rewrote function `PedCompare()` to increase clarity of code for easier maintenance; changed output format somewhat. 
- Added a vignette about the ageprior, and rewrote sections of the main vignette to incorporate new functions
- In the Fortran part, re-implemented how (candidate) (grand)parent-pairs are filtered and assigned 


### Bug fixes
- `GenoConvert()` skipped first individual when reading .raw file. Circumvent this bug in earlier versions by using option `header=FALSE` (then header row is removed only once...)
- `ConfProb()` expected input parameter `nSim` to be strictly integer, now relaxed to any value convertible to a whole number
- fixed `ERROR! ***Invalid ParProb!***` triggered when some SNPs are monomorphic
- fixed `SEGFAULT` triggered when some SNPs have very high missingness (>80%); possibly sibship size out of bounds
- fixed `Error arguments imply differing number of rows` when there are dummy parents of 1 sex only
- fixed various mostly minor bugs in the Fortran code
- fixed bugs regarding 'link time optimisation'



### Minor changes
- `LifeHistData` may have 2 additional columns, with minimum and maximum possible birth year
- second example pedigree (`Ped_griffin`) to illustrate overlapping generations, used in age vignette 
- `SummarySeq()`: added a pedigree summary table identical to a subset of the table returned by R package `pedantics`' `pedStatSummary`; that package has been archived on CRAN. Added option `Panels` to only plot (a) specific panel(s).  


# Sequoia 1.3.3

### Bug fixes
- fixes bug that caused R to crash (Fortran array indexing out-of-bounds)


# Sequoia 1.3.1

### New features
- several functions have become user-visible: `CheckGeno()`, `MkGenoErrors()`, `GetMaybeRel()`, `GetRelCat()`
- plotting functions added: `PlotAgePrior()` and `SummarySeq()`
- function `SimGeno()` input parameters have changed, old ones will be deprecated


### Minor changes
- extended vignette with function overview & FAQ
- numerous edits to fortran source code to better handle certain (rarer) types of relatives


### Bug fixes
- various bug fixes in fortran source code


# sequoia 1.1.1

### Major changes
- possibly. 


# sequoia 1.0.0

### New features
- added functions `EstConf`, `SnpStats`

### Major changes
- considerable changes in Fortran code


# sequoia 0.9.2

### New features
- added functionality for hermaphrodites (in silico cloned into male + female)


# sequoia 0.7.2
First version on CRAN!
