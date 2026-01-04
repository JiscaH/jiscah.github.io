# Convert Genotype Data

Convert genotype data in various formats to sequoia's
1-column-per-marker format, PLINK's ped format, or Colony's
2-columns-per-marker format.

## Usage

``` r
GenoConvert(
  InData = NULL,
  InFile = NULL,
  InFormat = "raw",
  OutFile = NA,
  OutFormat = "seq",
  Missing = c("-9", "NA", "??", "?", "NULL", "-1", c("0")[InFormat %in% c("col",
    "ped")]),
  sep = c(" ", "\t", ",", ";"),
  header = NA,
  IDcol = NA,
  FIDcol = NA,
  FIDsep = "__",
  dropcol = NA,
  quiet = FALSE
)
```

## Arguments

- InData:

  dataframe, matrix or
  [`genlight`](https://rdrr.io/pkg/adegenet/man/genlight.html) object
  with genotypes to be converted.

- InFile:

  character string with name of genotype file to be converted.

- InFormat:

  One of 'seq' (sequoia), 'ped' (PLINK .ped file), 'col' (COLONY), 'raw'
  (PLINK –recodeA), 'vcf' (requires library `{vcfR}`), 'single' (1
  column per SNP), or 'double' (2 columns per SNP); see Details.

- OutFile:

  character string with name of converted file. If NA, return matrix
  with genotypes in console (default); if NULL, write to
  'GenoForSequoia.txt' in current working directory.

- OutFormat:

  as `InFormat`; only 'seq', 'col', and 'ped' are implemented. For 'ped'
  also a sham .map file is created, so that the file can be read by
  PLINK. Only for 'ped' are extensions .ped & .map added to the
  specified OutFile filename.

- Missing:

  vector with symbols interpreted as missing data. '0' is missing data
  for `InFormat`s 'col' and 'ped' only.

- sep:

  vector with field separator strings that will be tried on `InFile`.
  Ignored if package data.table is present or if `InFormat='vcf'` or
  'vcf.gz'. The `OutFile` separator uses the
  [`write.table`](https://rdrr.io/r/utils/write.table.html) default,
  i.e. one blank space.

- header:

  a logical value indicating whether the file contains a header as its
  first line. If NA (default), set to TRUE for 'raw', and FALSE
  otherwise.

- IDcol:

  number giving the column with individual IDs; 0 indicates the rownames
  (for InData only). If NA (default), set to 2 for InFormat 'raw' and
  'ped', and otherwise to 1 for InFile and 0 (rownames) for InData,
  except when InData has a column labeled 'ID'.

- FIDcol:

  column with the family IDs, if any are wished to be used. This is
  column 1 for InFormat 'raw' and 'seq', but those are by default not
  used.

- FIDsep:

  string used to paste FID and IID together into a composite-ID (value
  passed to [`paste`](https://rdrr.io/r/base/paste.html)'s `collapse`).
  This joining can be reversed using
  [`PedStripFID`](https://jiscah.github.io/reference/PedStripFID.md).

- dropcol:

  columns to exclude from the output data, on top of IDcol and FIDcol
  (which become rownames). When NA, defaults to columns 3-6 for InFormat
  'raw' and 'seq'. Can also be used to drop some SNPs, see example below
  on how to do this for the 2-columns-per-SNP input formats.

- quiet:

  suppress messages and warnings.

## Value

A genotype matrix in the specified output format; the default sequoia
format ('seq') has 1 column per SNP coded in 0/1/2 format (major
homozygote /heterozygote /minor homozygote) with -9 for missing values,
sample IDs in row names and SNP names in column names. If 'OutFile' is
specified, the matrix is written to this file and nothing is returned
inside R.

## Details

The first two arguments are interchangeable, and can be given unnamed.
The first argument is assumed to be a file name if it is of class
'character' and length 1, and to be the genetic data if it is a matrix
or dataframe.

If package data.table is detected,
[`fread`](https://rdatatable.gitlab.io/data.table/reference/fread.html)
is used to read in the data from file. Otherwise, a combination of
[`readLines`](https://rdrr.io/r/base/readLines.html) and
[`strsplit`](https://rdrr.io/r/base/strsplit.html) is used.

## Input formats

The following formats can be specified by `InFormat`:

- seq:

  (sequoia) genotypes are coded as 0, 1, 2, missing as \\-9\\ (in input
  any negative number or NA are OK), in 1 column per marker. Column 1
  contains IDs, there is no header row.

- ped:

  (PLINK) genotypes are coded as A, C, T, G, missing as 0, in 2 columns
  per marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to 6
  ignored). If an associated .map file exists, SNP names will be read
  from there.

- raw:

  (PLINK) genotypes are coded as 0, 1, 2, missing as NA, in 1 column per
  marker. The first 6 columns are descriptive (1:FID, 2:IID, 3 to 6
  ignored), and there is a header row. This is produced by PLINK's
  option –recodeA

- col:

  (Colony) genotypes are coded as numeric values, missing as 0, in 2
  columns per marker. Column 1 contains IDs.

- vcf:

  (VCF) genotypes are coded as '0/0','0/1','1/1', variable number of
  header rows followed by 1 row per SNP, with various columns of
  metadata followed by 1 column per individual. Requires package vcfR.

- single:

  1 column per marker, otherwise unspecified

- double:

  2 columns per marker, otherwise unspecified

For each `InFormat`, its default values for
`Missing, header, IDcol, FIDcol`, and `dropcol` can be overruled by
specifying the corresponding input parameters.

## Error messages

Occasionally when reading in a file `GenoConvert` may give an error that
'rows have unequal length'. `GenoConvert` makes use of
[`readLines`](https://rdrr.io/r/base/readLines.html) and
[`strsplit`](https://rdrr.io/r/base/strsplit.html), which is much faster
than [`read.table`](https://rdrr.io/r/utils/read.table.html) for large
datafiles, but also more sensitive to unusual line endings, unusual
end-of-file characters, or invisible characters (spaces or tabs) after
the end of some lines. In these cases, try to read the data from file
using read.table or read.csv, and then use `GenoConvert` on this
dataframe or matrix, see example.

Any warnings generated by
[`CheckGeno`](https://jiscah.github.io/reference/CheckGeno.md) regarding
SNPs scored for few individuals and/or individuals scored for few SNPs
etc. are only for your information; none are excluded from
`GenoConvert`'s output, but these SNPs and/or individuals will be
excluded during pre-processing of the data in any of the other functions
in this package.

## See also

[`CheckGeno`](https://jiscah.github.io/reference/CheckGeno.md)`, `[`SnpStats`](https://jiscah.github.io/reference/SnpStats.md)`, `[`LHConvert`](https://jiscah.github.io/reference/LHConvert.md).

## Author

Jisca Huisman, <jisca.huisman@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires PLINK installed & in system PATH:

# tinker with window size, window overlap and VIF to get a set of
# 400 - 800 markers (100-200 enough for just parentage):
system("cmd", input = "plink --file mydata --indep 50 5 2")
system("cmd", input = "plink --file mydata --extract plink.prune.in
  --recodeA --out PlinkOUT")

GenoM <- GenoConvert(InFile = "PlinkOUT.raw", InFormat='raw')
# which is the same as
GenoM <- GenoConvert(PlinkOUT.raw, InFormat='single',
                    IDcol=2, dropcol=c(1,3:6), header=TRUE)
# (but it will complain that InFormat='single' is not consistent with .raw
# file extension)

# save time on file conversion next time:
write.table(GenoM, file="Geno_sequoia.txt", quote=FALSE, col.names=FALSE)
GenoM <- as.matrix(read.table("Geno_sequoia.txt", row.names=1, header=FALSE))

# drop some SNPs, e.g. after a warning of >2 alleles:
dropSNP <- c(5,68,101,128)
GenoM <- GenoConvert(ColonyFile, InFormat = "col",
                     dropcol = 1 + c(2*dropSNP-1, 2*dropSNP) )

# circumvent a 'rows have unequal length' error:
GenoTmp <- as.matrix(read.table("mydata.txt", header=TRUE, row.names=1))
GenoM <- GenoConvert(InData=GenoTmp, InFormat="single", IDcol=0)

# can also write to file, e.g. simulated genotypes:
GenoConvert(Geno_A, InFormat='seq', OutFormat='ped', OutFile = sim_genotypes)
} # }
```
