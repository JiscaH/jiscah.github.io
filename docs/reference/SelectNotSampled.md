# select non-genotyped parents

select non-genotyped parents

## Usage

``` r
SelectNotSampled(Ped, ParMis)
```

## Arguments

- Ped:

  pedigree, after PedPolish()

- ParMis:

  single number or vector length two with proportion of parents with
  fully missing genotype

## Value

vector with genotype matrix row numbers of non-sampled individuals
