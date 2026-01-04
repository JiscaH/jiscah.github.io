# Count Generations

For each individual in a pedigree, count the number of generations since
its most distant pedigree founder.

## Usage

``` r
getGenerations(Ped, StopIfInvalid = TRUE)
```

## Arguments

- Ped:

  dataframe, pedigree with the first three columns being id - dam -
  sire. Column names are ignored, as are additional columns.

- StopIfInvalid:

  if a pedigree loop is detected, stop with an error (TRUE, default) or
  return the Pedigree, to see where the problem(s) occur.

## Value

A vector with the generation number for each individual, starting at 0
for founders. Offspring of G0 X G0 are G1, offspring of G0 X G1 or G1 x
G1 are G2, etc. `NA` indicates a pedigree loop where an individual is
its own ancestor (or that the pedigree has \>1000 generations).

If no output name is specified, no results are returned, only an error
message when the pedigree contains a loop.

To get more details about a pedigree loop, you can use
https://github.com/JiscaH/sequoiaExtra/blob/main/find_pedigree_loop.R

## See also

[`GetAncestors`](https://jiscah.github.io/reference/GetAncestors.md)`, `[`GetDescendants`](https://jiscah.github.io/reference/GetDescendants.md)
to get all ancestors resp. descendants of a specific individual (with a
warning if it is its own ancestor);
[FindFamilies](https://jiscah.github.io/reference/FindFamilies.md) to
find connected sub-pedigrees.

## Examples

``` r
# returns nothing if OK, else error:
getGenerations(SeqOUT_griffin$Pedigree)

# returns vector with generation numbers:
G <- getGenerations(SeqOUT_griffin$Pedigree, StopIfInvalid=FALSE)
table(G, useNA='ifany')
#> G
#>  0  1  2  3  4  5  6  7  8 
#> 36 20 18 18 14 12 22 22  5 
Ped_plus_G <- cbind(SeqOUT_griffin$Pedigree, G)
```
