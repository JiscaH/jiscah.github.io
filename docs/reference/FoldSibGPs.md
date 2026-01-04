# Fold IDs of Sibship Grandparents

Fold IDs of sibship grandparents into a 2 x nInd/2 x 2 array, as they
are stored in Fortran, and then stretch this into a vector that can be
passed to Fortran and easily be transformed back into said 3D array.

## Usage

``` r
FoldSibGPs(PedNum, Ng, Nd)
```

## Arguments

- PedNum:

  pedigree, ids replaced by numbers, dummies negative.

- Ng:

  no. genotyped indivs.

- Nd:

  length 2 vector, no. female & male dummies.

## Value

An integer vector, with missing values as 0.
