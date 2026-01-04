# Plot Summary Overview of sequoia Output

visualise the numbers of assigned parents, sibship sizes, and parental
LLRs

## Usage

``` r
PlotSeqSum(SeqSum, Pedigree = NULL, Panels = "all", ask = TRUE)
```

## Arguments

- SeqSum:

  list output from
  [`SummarySeq`](https://jiscah.github.io/reference/SummarySeq.md).

- Pedigree:

  dataframe with at least id, dam and sire in columns 1-3, respectively.
  If columns with parental LLRs and/or Mendelian errors are present,
  these will be plotted as well.

- Panels:

  character vector with panel(s) to plot. Choose from 'all', 'G.parents'
  (parents of genotyped individuals), 'D.parents' (parents of dummies),
  'O.parents' (parents of non-genotyped non-dummies), sibships', 'LLR',
  'OH'.

- ask:

  ask for user key stroke before proceeding to next plot.

## Examples

``` r
sumry <- SummarySeq(SeqOUT_griffin, Plot=FALSE)
PlotSeqSum(sumry, SeqOUT_griffin$Pedigree, Panels='all', ask=FALSE)





```
