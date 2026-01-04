# Inheritance patterns

Inheritance patterns used by SimGeno for non-autosomal SNPs, identical
to those in Inherit.xlsx

## Usage

``` r
data(Inherit_patterns)
```

## Format

An array with the following dimensions:

- d1:

  type: autosomal, x-chromosome, y-chromosome, or mtDNA

- d2:

  offspring sex: female, male, or unknown

- d3:

  offspring genotype: aa (0), aA (1), Aa (1), or AA (2)

- d4:

  mother genotype

- d5:

  father genotype

## See also

[`SimGeno`](https://jiscah.github.io/reference/SimGeno.md)

## Author

Jisca Huisman, <jisca.huisman@gmail.com>
