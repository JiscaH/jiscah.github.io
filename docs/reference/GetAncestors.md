# Get ancestors

get all ancestors of an individual

## Usage

``` r
GetAncestors(id, Pedigree)
```

## Arguments

- id:

  id of the individual

- Pedigree:

  dataframe with columns id - parent1 - parent2; only the first 3
  columns will be used.

## Value

a list with as first element `id`, second parents, third grandparents,
etc.. Each element is a vector with ids, the first three elements are
named, the rest numbered. Ancestors are unsorted within each list
element.

## Examples

``` r
Anc_i200  <- GetAncestors('i200_2010_F', Ped_griffin)

```
