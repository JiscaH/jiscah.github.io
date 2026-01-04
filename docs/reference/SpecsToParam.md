# Specs to PARAM

Convert 1-row dataframe `Specs` into list `PARAM`, optionally including
various other objects in the list.

## Usage

``` r
SpecsToParam(Specs, ErrM = NULL, ErrFlavour = NULL, dimGeno = NULL, ...)
```

## Arguments

- Specs:

  1-row dataframe, element of
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md) output
  list.

- ...:

  other objects to append to the list, such as `ErrM` and `quiet`.

## Value

A named list.
