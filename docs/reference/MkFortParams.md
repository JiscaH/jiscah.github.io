# PARAM to FortPARAM

Convert list `PARAM` into a list with integer-only and double-only
vectors, to be passed to Fortran.

## Usage

``` r
MkFortParams(PARAM, fun = "main")
```

## Arguments

- PARAM:

  list with input parameters.

- fun:

  function from which `MkFortParams` is called, determines which
  elements are included in the output list.

## Value

A list with elements

- Ng:

  Integer, number of individuals

- SpecsInt:

  8 integers:

  - nSnp

  - MaxMismatchV; DUP - OH - ME

  - MaxSibshipSize

  - Complx, 0=mono, 1=simp, 2=full

  - quiet, -1=verbose, 0=FALSE, 1=TRUE

  - nAgeCl, nrow(AgePriors)

- SpecsDbl:

  2 double precision numbers:

  - Tfilter (\< 0)

  - Tassign (\> 0)

- ErrM:

  double, 3x3 matrix passed as length-9 vector

- SpecsIntMkPed:

  `fun='main'` only

  - AgeEffect, 0=no, 1=yes, 2=extra

  - CalcLLR, 0=FALSE, 1=TRUE

  - Herm, 0=no, 1= dam/sire distinction, 2=no dam/sire distinction

- SpecsIntAmb:

  `fun='mayberel'` only

  - ParSib 1=par, 2=ped

  - nAmbMax
