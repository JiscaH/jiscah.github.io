# Smooth out dips in ageprior matrix

...

## Usage

``` r
SmoothAP(V, tiny = 0.001)
```

## Arguments

- V:

  column in ageprior matrix (vector); strictly positive

- tiny:

  smallest non-zero value in V

## Details

Sets dips (\<10% of average of neighbouring ages) to the average of the
neighbouring ages, sets the age after the end (oldest observed age) to
LR(end)/2, and assigns a small value (0.001) to the ages before the
front (youngest observed age) and after the new end. Peaks are not
smoothed out, as these are less likely to cause problems than dips, and
are more likely to be genuine characteristics of the species.
