# Example output from check for relatives: griffins

Example output of a check for parent-offspring pairs and
parent-parent-offspring trios with
[`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md), with
[`Geno_griffin`](https://jiscah.github.io/reference/Geno_griffin.md) as
input (simulated from
[`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md)).

## Usage

``` r
data(MaybeRel_griffin)
```

## Format

a list with 2 dataframes, 'MaybePar' and 'MaybeTrio'. See
[`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md) for
further details.

## See also

[`SeqOUT_griffin`](https://jiscah.github.io/reference/SeqOUT_griffin.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
MaybeRel_griffin <- GetMaybeRel(GenoM = Geno_griffin, Err=0.001,
                                Module = 'par')
} # }
```
