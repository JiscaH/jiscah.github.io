# Example genotype file: 'HSg5'

Simulated genotype data for all\* individuals in Pedigree
[`Ped_HSg5`](https://jiscah.github.io/reference/Ped_HSg5.md) (\*: with
40

## Usage

``` r
data(Geno_HSg5)
```

## Format

A genotype matrix with 920 rows (ids) and 200 columns (SNPs). Each SNP
is coded as 0/1/2 copies of the reference allele, with -9 for missing
values. Ids are stored as rownames.

## See also

[`LH_HSg5`](https://jiscah.github.io/reference/LH_HSg5.md)`, `[`SimGeno`](https://jiscah.github.io/reference/SimGeno.md)`, `[`SeqOUT_HSg5`](https://jiscah.github.io/reference/SeqOUT_HSg5.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# this output was created as follows:
Geno_HSg5 <- SimGeno(Ped = Ped_HSg5, nSnp = 200, ParMis=0.4,
                     CallRate = 0.9, SnpError = 0.005)
} # }
```
