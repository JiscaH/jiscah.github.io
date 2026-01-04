# PARAM to Specs

Convert list `PARAM` into 1-row dataframe `Specs`. Only to be called by
[`sequoia`](https://jiscah.github.io/reference/sequoia.md).

## Usage

``` r
ParamToSpecs(PARAM, TimeStart, ErrFlavour)
```

## Arguments

- PARAM:

  list with input parameters.

- TimeStart:

  time at which
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md) run was
  started.

- ErrFlavour:

  character name or function.

## Value

The 1-row `Specs` dataframe.
