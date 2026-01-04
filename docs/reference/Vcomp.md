# Compare two vectors

Compare a vector with inferred sibs to a vector of \`true' sibs

## Usage

``` r
Vcomp(Infrd, Simld, SNPd)
```

## Arguments

- Infrd:

  vector of inferred sibs

- Simld:

  vector of true sibs

- SNPd:

  character vector with IDs of genotyped individuals

## Value

a named numeric vector of length 4, with the total length of Simld, the
length of the intersect of the two vectors, the number occurring in
Infrd but not Simld ('err'), and the number occuring in Simld but not
Infrd ('missed').
