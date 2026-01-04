# Generate Genotyping Error Matrix

Make a vector or matrix specifying the genotyping error pattern, or a
function to generate such a vector/matrix from a single value Err.

with the probabilities of observed genotypes (columns) conditional on
actual genotypes (rows), or return a function to generate such matrices
(using a single value Err as input to that function).

## Usage

``` r
ErrToM(Err = NA, flavour = "version2.9", Return = "matrix")
```

## Arguments

- Err:

  estimated genotyping error rate, as a single number, or 3x3 or 4x4
  matrix, or length 3 vector. If a single number, an error model is used
  that aims to deal with scoring errors typical for SNP arrays. If a
  matrix, this should be the probability of observed genotype (columns)
  conditional on actual genotype (rows). Each row must therefore sum
  to 1. If `Return='function'`, this may be `NA`. If a vector, these are
  the probabilities (observed given actual) hom\|other hom, het\|hom,
  and hom\|het.

- flavour:

  vector-generating or matrix-generating function, or one of
  'version2.9', 'version2.0', 'version1.3' (='SNPchip'), 'version1.1'
  (='version111'), referring to the sequoia version in which it was used
  as default. Only used if `Err` is a single number.

- Return:

  output, 'matrix' (default), 'vector', 'function' (matrix-generating),
  or 'v_function' (vector-generating)

## Value

Depending on `Return`, either:

- `'matrix'`: a 3x3 matrix, with probabilities of observed genotypes
  (columns) conditional on actual (rows)

- `'function'`: a function taking a single value `Err` as input, and
  generating a 3x3 matrix

- `'vector'`: a length 3 vector, with the probabilities (observed given
  actual) hom\|other hom, het\|hom, and hom\|het.

## Details

By default (`flavour` = "version2.9"), `Err` is interpreted as a
locus-level error rate (rather than allele-level), and equals the
probability that an actual heterozygote is observed as either homozygote
(i.e., the probability that it is observed as AA = probability that
observed as aa = `Err`/2). The probability that one homozygote is
observed as the other is (`Err`/2\\)^2\\.

The inbuilt 'flavours' correspond to the presumed and simulated error
structures, which have changed with sequoia versions. The most
appropriate error structure will depend on the genotyping platform;
'version0.9' and 'version1.1' were inspired by SNP array genotyping
while 'version1.3' and 'version2.0' are intended to be more general.

This function, and throughout the package, it is assumed that the two
alleles \\A\\ and \\a\\ are equivalent. Thus, using notation
\\P\\(observed genotype \|actual genotype), that \\P(AA\|aa) =
P(aa\|AA)\\, \\P(aa\|Aa)=P(AA\|Aa)\\, and \\P(aA\|aa)=P(aA\|AA)\\.

|             |              |               |              |
|-------------|--------------|---------------|--------------|
| **version** | **hom\|hom** | **het\|hom**  | **hom\|het** |
| **2.9**     | \\(E/2)^2\\  | \\E-(E/2)^2\\ | \\E/2\\      |
| **2.0**     | \\(E/2)^2\\  | \\E(1-E/2)\\  | \\E/2\\      |
| **1.3**     | \\(E/2)^2\\  | \\E\\         | \\E/2\\      |
| **1.1**     | \\E/2\\      | \\E/2\\       | \\E/2\\      |
| **0.9**     | \\0\\        | \\E\\         | \\E/2\\      |

or in matrix form, Pr(observed genotype (columns) \| actual genotype
(rows)):

*version2.9:*

|       |             |                |             |
|-------|-------------|----------------|-------------|
|       | **0**       | **1**          | **2**       |
| **0** | \\1-E\\     | \\E -(E/2)^2\\ | \\(E/2)^2\\ |
| **1** | \\E/2\\     | \\1-E\\        | \\E/2\\     |
| **2** | \\(E/2)^2\\ | \\E -(E/2)^2\\ | \\1-E\\     |

*version2.0:*

|       |               |              |               |
|-------|---------------|--------------|---------------|
|       | **0**         | **1**        | **2**         |
| **0** | \\(1-E/2)^2\\ | \\E(1-E/2)\\ | \\(E/2)^2\\   |
| **1** | \\E/2\\       | \\1-E\\      | \\E/2\\       |
| **2** | \\(E/2)^2\\   | \\E(1-E/2)\\ | \\(1-E/2)^2\\ |

*version1.3*

|       |                 |         |                 |
|-------|-----------------|---------|-----------------|
|       | **0**           | **1**   | **2**           |
| **0** | \\1-E-(E/2)^2\\ | \\E\\   | \\(E/2)^2\\     |
| **1** | \\E/2\\         | \\1-E\\ | \\E/2\\         |
| **2** | \\(E/2)^2\\     | \\E\\   | \\1-E-(E/2)^2\\ |

*version1.1*

|       |         |         |         |
|-------|---------|---------|---------|
|       | **0**   | **1**   | **2**   |
| **0** | \\1-E\\ | \\E/2\\ | \\E/2\\ |
| **1** | \\E/2\\ | \\1-E\\ | \\E/2\\ |
| **2** | \\E/2\\ | \\E/2\\ | \\1-E\\ |

*version0.9* (not recommended)

|       |         |         |         |
|-------|---------|---------|---------|
|       | **0**   | **1**   | **2**   |
| **0** | \\1-E\\ | \\E\\   | \\0\\   |
| **1** | \\E/2\\ | \\1-E\\ | \\E/2\\ |
| **2** | \\0\\   | \\E\\   | \\1-E\\ |

When `Err` is a length 3 vector, or if `Return = 'vector'` these are the
following probabilities:

- hom\|hom: an actual homozygote is observed as the other homozygote
  (\\E_1\\)

- het\|hom: an actual homozygote is observed as heterozygote (\\E_2\\)

- hom\|het: an actual heterozygote is observed as homozygote (\\E_3\\)

and Pr(observed genotype (columns) \| actual genotype (rows)) is then:

|       |               |            |               |
|-------|---------------|------------|---------------|
|       | **0**         | **1**      | **2**         |
| **0** | \\1-E_1-E_2\\ | \\E_2\\    | \\E_1\\       |
| **1** | \\E_3\\       | \\1-2E_3\\ | \\E_3\\       |
| **2** | \\E_1\\       | \\E_2\\    | \\1-E_1-E_2\\ |

When the SNPs are scored via sequencing (e.g. RADseq or DArTseq), the
3rd error rate (hom\|het) is typically considerably higher than the
other two, while for SNP arrays it tends to be similar to P(het\|hom).

## Examples

``` r
ErM <- ErrToM(Err = 0.05)
ErM
#>    obs
#> act        0        1        2
#>   0 0.950000 0.049375 0.000625
#>   1 0.025000 0.950000 0.025000
#>   2 0.000625 0.049375 0.950000
ErrToM(ErM, Return = 'vector')
#>  hom|hom  het|hom  hom|het 
#> 0.000625 0.049375 0.025000 


# use error matrix from Whalen, Gorjanc & Hickey 2018
funE <- function(E) {
 matrix(c(1-E*3/4, E/2, E/4,
          E/4, 1-2*E/4, E/4,
          E/4, E/2, 1-E*3/4),
          3,3, byrow=TRUE)  }
ErrToM(Err = 0.05, flavour = funE)
#>        [,1]  [,2]   [,3]
#> [1,] 0.9625 0.025 0.0125
#> [2,] 0.0125 0.975 0.0125
#> [3,] 0.0125 0.025 0.9625
# equivalent to:
ErrToM(Err = c(0.05/4, 0.05/2, 0.05/4))
#>    obs
#> act      0     1      2
#>   0 0.9625 0.025 0.0125
#>   1 0.0125 0.975 0.0125
#>   2 0.0125 0.025 0.9625
```
