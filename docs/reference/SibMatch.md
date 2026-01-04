# Find the closest matching inferred sibship to a true sibship

Find the closest matching inferred sibship to a true sibship

## Usage

``` r
SibMatch(SimX, Infrd, SNPd)
```

## Arguments

- SimX:

  a vector with the IDs in the true (Ped1) sibship

- Infrd:

  a list of vectors with the IDs in the inferred (Ped2) sibships

- SNPd:

  character vector with IDs of genotyped individuals

## Value

a named numeric vector with the number of matches ('NumMatch'), the
position of the best match ('Best'), the inferred sibship size of this
best match ('Tot'), the number of matching IDs ('OK'), and the number of
mismatches ('err').
