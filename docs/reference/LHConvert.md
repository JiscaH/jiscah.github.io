# Extract Sex and Birth Year from PLINK File

Convert the first six columns of a PLINK .fam, .ped or .raw file into a
three-column lifehistory file for sequoia. Optionally FID and IID are
combined.

## Usage

``` r
LHConvert(
  PlinkFile = NULL,
  UseFID = FALSE,
  SwapSex = TRUE,
  FIDsep = "__",
  LifeHistData = NULL
)
```

## Arguments

- PlinkFile:

  character string with name of genotype file to be converted.

- UseFID:

  use the family ID column. The resulting ids (rownames of GenoM) will
  be in the form FID\_\_IID.

- SwapSex:

  change the coding from PLINK default (1=male, 2=female) to sequoia
  default (1=female, 2=male); any other numbers are set to NA.

- FIDsep:

  characters inbetween FID and IID in composite-ID. By default a double
  underscore is used, to avoid problems when some IIDs contain an
  underscore. Only used when UseFID=TRUE.

- LifeHistData:

  dataframe with additional sex and birth year info. In case of
  conflicts, LifeHistData takes priority, with a warning. If
  UseFID=TRUE, IDs in LifeHistData are assumed to be already as
  FID\_\_IID.

## Value

A dataframe with id, sex and birth year, which can be used as input for
[`sequoia`](https://jiscah.github.io/reference/sequoia.md).

## Details

The first 6 columns of PLINK .fam, .ped and .raw files are by default
FID - IID - father ID (ignored) - mother ID (ignored) - sex - phenotype.

## See also

[`GenoConvert`](https://jiscah.github.io/reference/GenoConvert.md),
[`PedStripFID`](https://jiscah.github.io/reference/PedStripFID.md) to
reverse `UseFID`.

## Examples

``` r
if (FALSE) { # \dontrun{
# combine FID and IID in dataframe with additional sex & birth years
ExtraLH$FID_IID <- paste(ExtraLH$FID, ExtraLH$IID, sep = "__")
LH.new <- LHConvert(PlinkFile, UseFID = TRUE, FIDsep = "__",
                    LifeHistData = ExtraLH)
} # }
```
