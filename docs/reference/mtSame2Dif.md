# Check and recode mtSame matrix

Recode to 1=different, 0=same

## Usage

``` r
mtSame2Dif(mtSame = NULL, gID = NULL)
```

## Arguments

- mtSame:

  matrix indicating whether individuals (might) have the same
  mitochondrial haplotype (1), or definitely not (0). Not all
  individuals need to be included and order is not important, may not be
  square.

- gID:

  rownames of \`GenoM\`

## Value

square gID x gID matrix indicating whether individuals have definitely a
different mitochondrial haplotype (1), or (possibly) the same (0).
