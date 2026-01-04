# Find siblings

Find siblings

## Usage

``` r
rc(x, Ped)
```

## Arguments

- x:

  an ID

- Ped:

  a pedigree with columns id - dam - sire

## Value

The individuals which are full or half siblings to x, as a three-column
matrix with column names id1 (x), id2 (the siblings), and RC (the
relatedness category, 'FS' or 'HS').
