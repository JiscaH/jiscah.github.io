# Plot Age Priors

Visualise the age-difference based prior probability ratios as a
heatmap.

## Usage

``` r
PlotAgePrior(AP = NULL, legend = TRUE)
```

## Arguments

- AP:

  matrix with age priors (\\P(A\|R)/P(A)\\) with age differences in rows
  and relationships in columns; by default M: maternal parent (mother),
  P: paternal parent (father), FS: full siblings, MS: maternal siblings
  (full + half), PS: paternal siblings.

- legend:

  if `TRUE`, a new plotting window is started and
  [`layout`](https://rdrr.io/r/graphics/layout.html) is used to plot a
  legend next to the main plot. Set to `FALSE` if you want to add it as
  panel to an existing plot (e.g. with `par(mfcol=c(2,2))`).

## Value

A heatmap.

## See also

[`MakeAgePrior`](https://jiscah.github.io/reference/MakeAgePrior.md),
[`SummarySeq`](https://jiscah.github.io/reference/SummarySeq.md).

## Examples

``` r
PlotAgePrior(SeqOUT_griffin$AgePriors)

PlotAgePrior(SeqOUT_griffin$AgePriorExtra)

```
