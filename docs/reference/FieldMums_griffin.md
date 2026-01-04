# Example field-observed mothers: griffins

Example field pedigree used in vignette for
[`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)
example. Non-genotyped females have IDs 'BlueRed', 'YellowPink', etc.

## Usage

``` r
data(FieldMums_griffin)
```

## Format

A data frame with 144 rows and 2 variables (id, mum)

## See also

[`SeqOUT_griffin`](https://jiscah.github.io/reference/SeqOUT_griffin.md)
for a sequoia run on simulated genotype data,
[`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md) for
the 'true' pedigree.

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
PC_griffin <- PedCompare(Ped1 = cbind(FieldMums_griffin, sire=NA),
                         Ped2 = SeqOUT_griffin$Pedigree)
} # }
```
