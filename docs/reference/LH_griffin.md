# Example life history data: griffins

Example life history data associated with the griffin pedigree.

## Usage

``` r
data(LH_griffin)
```

## Format

A data frame with 200 rows and 3 variables (ID, Sex, BirthYear)

## See also

[`Ped_griffin`](https://jiscah.github.io/reference/Ped_griffin.md),
[`SeqOUT_griffin`](https://jiscah.github.io/reference/SeqOUT_griffin.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
BY <- rep(c(2001:2010), each=20)
Sex <- sample.int(n=2, size=200, replace=TRUE)
ID <- paste0("i", formatC(1:200, width=3, flag="0"), "_", BY, "_",
             ifelse(Sex==1, "F", "M"))
LH_griffin <- data.frame(ID, Sex, BirthYear = BY)
} # }
```
