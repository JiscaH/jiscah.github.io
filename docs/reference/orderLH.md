# Order Lifehistory Data

Order lifehistory data to match order of IDs in genotype data, filling
in gaps with missing values.

## Usage

``` r
orderLH(LH = NULL, gID = NULL)
```

## Arguments

- LH:

  dataframe with lifehistory information:

  ID

  :   max. 30 characters long,

  Sex

  :   1 = females, 2 = males, other numbers = unknown,

  Birth Year

  :   (or hatching year) Use negative numbers to denote missing values.

  BY.min

  :   minimum birth year (optional)

  BY.max

  :   maximum birth year (optional)

- gID:

  character vector with IDs in genotype data, in order of occurrence.

## Value

A dataframe with the same 5 columns, but with individuals in exactly the
same order as gID, including padding with 'empty' rows if an individual
in gID was not in the input-LH. Missing values are recoded to 3 for the
'Sex' column, and -999 for the birth year columns.
