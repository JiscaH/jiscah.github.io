# Special Merge

As regular merge, but combine data from columns with the same name.

## Usage

``` r
MergeFill(df1, df2, by, overwrite = FALSE, ...)
```

## Arguments

- df1:

  first dataframe (lowest priority if `overwrite=TRUE`).

- df2:

  second dataframe (highest priority if `overwrite=TRUE`).

- by:

  columns used for merging, required.

- overwrite:

  If FALSE (the default), NA's in df1 are replaced by values from df2.
  If TRUE, all values in df1 are overwritten by values from df2, except
  where df2 has NA.

- ...:

  additional arguments to merge, such as `all`.
