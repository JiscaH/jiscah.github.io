# Example output from pedigree inference: griffins

Example output of a sequoia run including sibship clustering, with
[`Geno_griffin`](https://jiscah.github.io/reference/Geno_griffin.md) as
input (simulated from
[`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md)).

## Usage

``` r
data(SeqOUT_griffin)
```

## Format

a list, see [`sequoia`](https://jiscah.github.io/reference/sequoia.md)

## See also

[`sequoia`](https://jiscah.github.io/reference/sequoia.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
SeqOUT_griffin <- sequoia(GenoM = Geno_griffin,
                          LifeHistData = LH_griffin,
                          Module = 'ped')
} # }
```
