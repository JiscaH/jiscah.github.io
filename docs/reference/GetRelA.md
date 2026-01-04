# Array with Pairwise Relationships

Generate an array indicating the relationship(s) between all pairs of
individuals according to the pedigree.

## Usage

``` r
GetRelA(Ped = NULL, GenBack = 1, patmat = TRUE, directed = TRUE, List = FALSE)
```

## Arguments

- Ped:

  dataframe with columns id - dam - sire.

- GenBack:

  number of generations back to consider; 1 returns parent-offspring and
  sibling relationships, 2 also returns grand-parental, avuncular and
  first cousins.

- patmat:

  logical, distinguish between paternal versus maternal relative pairs?
  For avuncular pairs, the distinction is never made.

- directed:

  logical, distinguish between 'O' vs 'P' or group into 'PO' ?

- List:

  logical, return a list instead of the default array

## Value

a 3D array indicating if the pair has the specified relationship (1) or
not (0). The various relationship considered are in the 3rd dimension:

- M:

- P:

- FS:

  full siblings, including double 'other half sibs'

- MS:

- PS:

- XS:

  other sibs: mother of A is father of B, or vv

etc.
