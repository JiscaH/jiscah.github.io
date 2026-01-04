# Write Sequoia Output to File

The various list elements returned by `sequoia` are each written to text
files in the specified folder, or to separate sheets in a single excel
file (requires library openxlsx).

## Usage

``` r
writeSeq(
  SeqList,
  GenoM = NULL,
  MaybeRel = NULL,
  PedComp = NULL,
  OutFormat = "txt",
  folder = "Sequoia-OUT",
  file = "Sequoia-OUT.xlsx",
  quiet = FALSE
)
```

## Arguments

- SeqList:

  list returned by
  [`sequoia`](https://jiscah.github.io/reference/sequoia.md), to be
  written out.

- GenoM:

  matrix with genetic data (optional). Ignored if OutFormat='xls', as
  the resulting file could become too large for excel.

- MaybeRel:

  list with results from
  [`GetMaybeRel`](https://jiscah.github.io/reference/GetMaybeRel.md)
  (optional).

- PedComp:

  list with results from
  [`PedCompare`](https://jiscah.github.io/reference/PedCompare.md)
  (optional). `SeqList$DummyIDs` is combined with `PedComp$DummyMatch`
  if both are provided.

- OutFormat:

  'xls' or 'txt'.

- folder:

  the directory where the text files will be written; will be created if
  it does not already exists. Relative to the current working directory,
  or NULL for current working directory. Ignored if `OutFormat='xls'`.

- file:

  the name of the excel file to write to, ignored if `OutFormat='txt'`.

- quiet:

  suppress messages.

## Details

The text files can be used as input for the stand-alone Fortran version
of sequoia, e.g. when the genotype data is too large for R. See
`vignette('sequoia')` for further details.

## See also

[`writeColumns`](https://jiscah.github.io/reference/writeColumns.md) to
write to a text file, using white space padding to keep columns aligned.

## Examples

``` r
if (FALSE) { # \dontrun{
writeSeq(SeqList, OutFormat="xls", file="MyFile.xlsx")

# add additional sheet to the excel file:
library(openxlsx)
wb <- loadWorkbook("MyFile.xlsx")
addWorksheet(wb, sheetName = "ExtraData")
writeData(wb, sheet = "ExtraData", MyData, rowNames=FALSE)
saveWorkbook(wb, "MyFile.xlsx", overwrite=TRUE, returnValue=TRUE)

# or: (package requires java & is trickier to install)
xlsx::write.xlsx(MyData, file = "MyFile.xlsx", sheetName="ExtraData",
      col.names=TRUE, row.names=FALSE, append=TRUE, showNA=FALSE)
} # }
```
