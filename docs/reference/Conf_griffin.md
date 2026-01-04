# Example output from estimating confidence probabilities: griffins

Example output of
[`EstConf`](https://jiscah.github.io/reference/EstConf.md), with the
inferred pedigree in
[`SeqOUT_griffin`](https://jiscah.github.io/reference/SeqOUT_griffin.md)
used as reference pedigree.

## Usage

``` r
data(Conf_griffin)
```

## Format

a list, see [`sequoia`](https://jiscah.github.io/reference/sequoia.md)

## See also

[`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md),
[`Geno_griffin`](https://jiscah.github.io/reference/Geno_griffin.md),

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
Conf_griffin <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
                        LifeHistData = LH_griffin,
                        args.sim = list(nSnp = 400, SnpError = 0.001,
                                        ParMis=0.4),
                        args.seq = list(Module = 'ped', Err=0.001),
                        nSim = 20,
                        nCores = 5,
                        quiet = TRUE)
} # }
```
