# Fix Pedigree

Ensure all parents & all genotyped individuals are included, remove
duplicates, rename columns, and replace 0 by NA or v.v.. Does not sort
parents before offspring.

## Usage

``` r
PedPolish(
  Pedigree,
  gID = NULL,
  ZeroToNA = TRUE,
  NAToZero = FALSE,
  DropNonSNPd = TRUE,
  addParentRows = TRUE,
  FillParents = FALSE,
  KeepAllColumns = TRUE,
  KeepAllRows = FALSE,
  NullOK = FALSE,
  LoopCheck = TRUE,
  StopIfInvalid = TRUE
)
```

## Arguments

- Pedigree:

  dataframe where the first 3 columns are id, dam, sire.

- gID:

  character vector with ids of genotyped individuals (rownames of
  genotype matrix).

- ZeroToNA:

  logical, replace 0's for missing values by NA's (defaults to `TRUE`).

- NAToZero:

  logical, replace NA's for missing values by 0's. If `TRUE`, ZeroToNA
  is automatically set to `FALSE`.

- DropNonSNPd:

  logical, remove any non-genotyped individuals (but keep non-genotyped
  parents), & sort pedigree in order of `gID`.

- addParentRows:

  add rows for any dams, sires, or individuals in `gID` not yet
  occurring in the id column.

- FillParents:

  logical, for individuals with only 1 parent assigned, set the other
  parent to a dummy (without assigning siblings or grandparents). Makes
  the pedigree compatible with R packages and software that requires
  individuals to have either 2 or 0 parents, such as
  [`kinship`](https://rdrr.io/pkg/kinship2/man/kinship.html).

- KeepAllColumns:

  Keep all columns in `Pedigree` (TRUE, default), or only id - dam -
  sire (FALSE).

- KeepAllRows:

  Keep all rows in `Pedigree` (TRUE), or drop rows where id = `NA`
  (FALSE, default). Duplicated rows are always removed.

- NullOK:

  logical, is it OK for Ped to be NULL? Then NULL will be returned.

- LoopCheck:

  logical, check for invalid pedigree loops by calling
  [`getGenerations`](https://jiscah.github.io/reference/getGenerations.md).

- StopIfInvalid:

  if a pedigree loop is detected, stop with an error (TRUE, default).

## Details

Recognized column names are an exact or partial match with (case is
ignored):

- id:

  "id", "iid", "off"

- dam:

  "dam", "mother", "mot", "mom", "mum", "mat"

- sire:

  "sire", "father", "fat", "dad", "pat"

`sequoia` requires the column order id - dam - sire; columns 2 and 3 are
swapped by this function if necessary.

## Examples

``` r
PedZ <- Ped_HSg5[41:50, ]
PedPolish(PedZ)
#>        id    dam   sire
#> 1  a00008   <NA>   <NA>
#> 2  a00013   <NA>   <NA>
#> 3  a00014   <NA>   <NA>
#> 4  a01001 a00014 b00011
#> 5  a01002 a00014 b00011
#> 6  a01005 a00013 b00001
#> 7  a01008 a00013 b00001
#> 8  a01010 a00008 b00016
#> 9  b00001   <NA>   <NA>
#> 10 b00011   <NA>   <NA>
#> 11 b00016   <NA>   <NA>
#> 12 b01003 a00014 b00011
#> 13 b01004 a00014 b00011
#> 14 b01006 a00013 b00001
#> 15 b01007 a00013 b00001
#> 16 b01009 a00008 b00016
PedPolish(PedZ, gID = rownames(SimGeno_example)[30:40], DropNonSNPd=TRUE)
#>        id    dam   sire
#> 11 b00009   <NA>   <NA>
#> 12 b00010   <NA>   <NA>
#> 10 b00007   <NA>   <NA>
#> 4  a01001 a00014 b00011
#> 5  a01002 a00014 b00011
#> 15 b01003 a00014 b00011
#> 16 b01004 a00014 b00011
#> 6  a01005 a00013 b00001
#> 17 b01006 a00013 b00001
#> 18 b01007 a00013 b00001
#> 7  a01008 a00013 b00001

if (FALSE) { # \dontrun{
# To get the output pedigree into kinship2 compatible format:
PedP <- sequoia::PedPolish(SeqOUT$Pedigree, DropNonSNPd=FALSE,
                           FillParents = TRUE)
PedP$Sex <- with(PedP, ifelse(id %in% dam, "female",  "male"))
# default to 'male' to avoid warning: "More than 25% of the gender values are
#  'unknown'"

Ped.fix <- with(PedP, kinship2::fixParents(id=id, dadid=sire, momid=dam,
                                           sex=Sex))
Ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))
} # }
```
