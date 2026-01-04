# Compare Pairwise Relationships

Compare, count and identify different types of relative pairs between
two pedigrees, or within one pedigree.

## Usage

``` r
ComparePairs(
  Ped1 = NULL,
  Ped2 = NULL,
  Pairs2 = NULL,
  GenBack = 1,
  patmat = FALSE,
  ExcludeDummies = TRUE,
  DumPrefix = c("F0", "M0"),
  Return = "Counts",
  Pairs_suffix = "?"
)
```

## Arguments

- Ped1:

  first (e.g. original/reference) pedigree, dataframe with 3 columns:
  id-dam-sire.

- Ped2:

  optional second (e.g. inferred) pedigree.

- Pairs2:

  optional dataframe with as first three columns: ID1-ID2- relationship,
  e.g. as returned by
  [`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md).
  Column names and any additional columns are ignored. May be provided
  in addition to, or instead of `Ped2`.

- GenBack:

  number of generations back to consider; 1 returns parent-offspring and
  sibling relationships, 2 also returns grandparental, avuncular and
  first cousins. GenBack \>2 is not implemented.

- patmat:

  logical, distinguish between paternal versus maternal relative pairs?

- ExcludeDummies:

  logical, exclude dummy IDs from output? Individuals with e.g. the same
  dummy father will still be counted as paternal halfsibs. No attempt is
  made to match dummies in one pedigree to individuals in the other
  pedigree; for that use
  [`PedCompare`](https://jiscah.github.io/reference/PedCompare.md).

- DumPrefix:

  character vector with the prefixes identifying dummy individuals. Use
  'F0' ('M0') to avoid matching to regular individuals with IDs starting
  with 'F' ('M'), provided `Ped2` has fewer than 999 dummy females
  (males).

- Return:

  return a matrix with `Counts` or a `Summary` of the number of
  identical relationships and mismatches per relationship, or detailed
  results as a 2xNxN `Array` or as a `Dataframe`. `All` returns a list
  with all four.

- Pairs_suffix:

  symbol added to the relationship abbreviations derived from `Pairs2`,
  when both `Ped2` and `Pairs2` are provided. Can be an empty string.

## Value

Depending on `Return`, one of the following, or a list with all:

- Counts:

  (the default), a matrix with counts, with the classification in `Ped1`
  on rows and that in `Ped2` in columns. Counts for 'symmetrical' pairs
  ("FS", "HS", "MHS", "PHS", "FC1", "DFC1", "U","X") are divided by two.

- Summary:

  a matrix with one row per relationship type and four columns , named
  as if `Ped1` is the true pedigree:

  n

  :   total number of pairs with that relationship in `Ped1`, and
      occurring in `Ped2`

  OK

  :   Number of pairs with same relationship in `Ped2` as in `Ped1`

  hi

  :   Number of pairs with 'higher' relationship in `Ped2` as in `Ped1`
      (e.g. FS instead of HS; ranking is the order given below)

  lo

  :   Number of pairs with 'lower' relationship in `Ped2` as in `Ped1`,
      but not unrelated in `Ped2`

- Array:

  a 2xNxN array (if `Ped2` or `Pairs2` is specified) or a NxN matrix ,
  where N is the total number of individuals occurring in `Ped1` and/or
  `Ped2`.

- Dataframe:

  a dataframe with \\N^2\\ rows and four columns:

  id.A

  :   First individual of the pair

  id.B

  :   Second individual of the pair

  RC1

  :   the relationship category in `Ped1`, as a factor with all
      considered categories as levels, including those with 0 count

  RC2

  :   the relationship category in `Ped2`

  Each pair is listed twice, e.g. once as P and once as O, or twice as
  FS.

## Details

If `Pairs2` is as returned by
[`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md)
(identified by the additional column names 'LLR' and 'OH'), these
relationship categories are appended with an '?' in the output, to
distinguish them from those derived from `Ped2`.

When `Pairs2$TopRel` contains values other than the ones listed among
the return values for the combination of `patmat` and `GenBack`, they
are prioritised in decreasing order of factor levels, or in decreasing
alphabetical order, and before the default (`ped2` derived) levels.

The matrix returned by `DyadCompare` \[Deprecated\] is a subset of the
matrix returned here using default settings.

## Relationship abbreviations and ranking

By default (`GenBack=1, patmat=FALSE`) the following 7 relationships are
distinguished:

- **S**: Self (not included in `Counts`)

- **MP**: Parent

- **O**: Offspring (not included in `Counts`)

- **FS**: Full sibling

- **HS**: Half sibling

- **U**: Unrelated, or otherwise related

- **X**: Either or both individuals not occurring in both pedigrees

In the array and dataframe, 'MP' indicates that the second (column)
individual is the parent of the first (row) individual, and 'O'
indicates the reverse.

When `GenBack=1, patmat=TRUE` the categories are (S)-M-P-(O)-FS-MHS-PHS-
U-X.

When `GenBack=2, patmat=TRUE`, the following relationships are
distinguished:

- **S**: Self (not included in `Counts`)

- **M**: Mother

- **P**: Father

- **O**: Offspring (not included in `Counts`)

- **FS**: Full sibling

- **MHS**: Maternal half-sibling

- **PHS**: Paternal half-sibling

- **MGM**: Maternal grandmother

- **MGF**: Maternal grandfather

- **PGM**: Paternal grandmother

- **PGF**: Paternal grandfather

- **GO**: Grand-offspring (not included in `Counts`)

- **FA**: Full avuncular; maternal or paternal aunt or uncle

- **HA**: Half avuncular

- **FN**: Full nephew/niece (not included in `Counts`)

- **HN**: Half nephew/niece (not included in `Counts`)

- **FC1**: Full first cousin

- **DFC1**: Double full first cousin

- **U**: Unrelated, or otherwise related

- **X**: Either or both individuals not occurring in both pedigrees

Note that for avuncular and cousin relationships no distinction is made
between paternal versus maternal, as this may differ between the two
individuals and would generate a large number of sub-classes. When a
pair is related via multiple paths, the first-listed relationship is
returned. To get all the different paths between a pair, use
[`GetRelM`](https://jiscah.github.io/reference/GetRelM.md) with
`Return='Array'`.

When `GenBack=2, patmat=FALSE`, MGM, MGF, PGM and PGF are combined into
GP, with the rest of the categories analogous to the above.

## See also

[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md) for
individual-based comparison;
[`GetRelM`](https://jiscah.github.io/reference/GetRelM.md) for a
pairwise relationships matrix of a single pedigree;
[`PlotRelPairs`](https://jiscah.github.io/reference/PlotRelPairs.md) for
visualisation of relationships within each pedigree.

To estimate P(actual relationship (Ped1) \| inferred relationship
(Ped2)), see examples at
[`EstConf`](https://jiscah.github.io/reference/EstConf.md).

## Examples

``` r
PairsG <- ComparePairs(Ped_griffin, SeqOUT_griffin[["Pedigree"]],
                       patmat = TRUE, ExcludeDummies = TRUE, Return = "All")
PairsG$Counts
#>      Ped2
#> Ped1     M    P    O   FS  MHS  PHS    U    X
#>   M     65    0    0    0    0    0    0  102
#>   P      0   79    0    0    0    0    0   84
#>   FS     0    0    0    5    0    0    0    0
#>   MHS    0    0    0    0   89    0    6  116
#>   PHS    0    0    0    0    0   76    6   72
#>   U      0    0    0    0    0    0 9685 9515
#>   X      0    0    0    0    0    0    0    0

# pairwise correct assignment rate:
PairsG$Summary[,"OK"] / PairsG$Summary[,"n"]
#>         M         P        FS       MHS       PHS         U 
#> 1.0000000 1.0000000 1.0000000 0.9368421 0.9268293 1.0000000 

# check specific pair:
PairsG$Array[, "i190_2010_M", "i168_2009_F"]
#> Ped1 Ped2 
#>  "M"  "X" 
# or
RelDF <- PairsG$Dataframe   # for brevity
RelDF[RelDF$id.A=="i190_2010_M" & RelDF$id.B=="i168_2009_F", ]
#>              id.A        id.B Ped1 Ped2
#> 33590 i190_2010_M i168_2009_F    M    X

# Colony-style lists of full sib dyads & half sib dyads:
FullSibDyads <- with(RelDF, RelDF[Ped1 == "FS" & id.A < id.B, ])
HalfSibDyads <- with(RelDF, RelDF[Ped1 == "HS" & id.A < id.B, ])
# Use 'id.A < id.B' because each pair is listed 2x
```
