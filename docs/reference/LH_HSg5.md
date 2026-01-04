# Example life history file: 'HSg5'

This is the life history file associated with
[`Ped_HSg5`](https://jiscah.github.io/reference/Ped_HSg5.md), which is
**Pedigree II** in the paper.

## Usage

``` r
data(LH_HSg5)
```

## Format

A data frame with 1000 rows and 3 variables:

- ID:

  Female IDs start with 'a', males with 'b'; the next 2 numbers give the
  generation number (00 – 05), the last 3 numbers the individual ID
  number (runs continuously across all generations)

- Sex:

  1 = female, 2 = male

- BirthYear:

  from 2000 (generation 0, founders) to 2005

## References

Huisman, J. (2017) Pedigree reconstruction from SNP data: Parentage
assignment, sibship clustering, and beyond. Molecular Ecology Resources
17:1009–1024.

## See also

[`Ped_HSg5`](https://jiscah.github.io/reference/Ped_HSg5.md)` `[`sequoia`](https://jiscah.github.io/reference/sequoia.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>
