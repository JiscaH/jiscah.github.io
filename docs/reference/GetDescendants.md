# Get descendants

get all descendants of an individual

## Usage

``` r
GetDescendants(id, Pedigree)
```

## Arguments

- id:

  id of the individual

- Pedigree:

  dataframe with columns id - parent1 - parent2; only the first 3
  columns will be used.

## Value

a list with as first element `id`, second offspring, third
grand-offspring, etc.. Each element is a vector with ids, the first
three elements are named, the rest numbered.
