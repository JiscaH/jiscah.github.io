# Scatter Plot of Pair LLRs

Plot LLR(rely/U) against LLR(relx/U), for one combination of
relationships, colour coded by fcl & top.

## Usage

``` r
LLRplot(relx, rely, LLRU, fcl, top, RelCol, bgcol, Tassign = 0.5, Tfilter = -2)
```

## Arguments

- relx:

  relationship to plot on the x-axis. One of 'PO', 'FS', 'HS', 'GP',
  'FA', or 'HA'.

- rely:

  relationship to plot on the y-axis; as `relx`.

- LLRU:

  matrix with log10-likelihoods, already scaled by LL(U) for each pair.

- fcl:

  focal relationship, sets outer circle colour of points.

- top:

  most likely relationship, sets inner filling colour of points.

- RelCol:

  named character vector with colours to use per relationship.

- bgcol:

  do background colour TRUE/FALSE.

- Tassign:

  assignment threshold, shown as grey square in bottom-left corner and
  band along the diagonal.

- Tfilter:

  filter threshold, shown as dark grey square in bottom-left.

## Details

The background of the plot is coloured to match `relx` (bottom triangle)
and `rely` (upper triangle).
