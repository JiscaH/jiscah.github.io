# Write Data to a File Column-wise

Write data.frame or matrix to a text file, using white space padding to
keep columns aligned as in `print`.

## Usage

``` r
writeColumns(x, file = "", row.names = TRUE, col.names = TRUE)
```

## Arguments

- x:

  the object to be written, preferably a matrix or data frame. If not,
  it is attempted to coerce x to a matrix.

- file:

  a character string naming a file.

- row.names:

  a logical value indicating whether the row names of x are to be
  written along with x.

- col.names:

  a logical value indicating whether the column names of x are to be
  written along with x.
