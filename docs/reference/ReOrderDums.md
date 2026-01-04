# Re-order likelihood vector

change the order for dummy individuals so that the results are easily
comparable with genotyped individuals

## Usage

``` r
ReOrderDums(LLM, Complex = "full")
```

## Arguments

- LLM:

  a matrix with log10-likelihoods; individuals in rows.

## Value

a matrix with similar dimensions, but changed column order:

- .:

  GP \\\rightarrow\\ PO (GP of sibship = parent of dummy)

- .:

  FA \\\rightarrow\\ FS

- .:

  HA \\\rightarrow\\ HS

- .:

  FS \\\rightarrow\\ FA (FS of sibs=offspring of dummy, NOT parent)

- .:

  HS \\\rightarrow\\ HA

- .:

  U \\\rightarrow\\ U

- .:

  PO \\\rightarrow\\ DUP (only if 'DUP' already among column names)
