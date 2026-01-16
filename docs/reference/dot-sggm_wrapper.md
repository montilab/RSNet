# A Wrapper function of SILGGM

A Wrapper function of SILGGM

## Usage

``` r
.sggm_wrapper(mat, method, alpha = 0.05, global = TRUE)
```

## Arguments

- mat:

  a n x p matrix, n = number of samples, and p = number of features

- method:

  the inference method

- alpha:

  A user-supplied sequence of pre-sepecified alpha levels for FDR
  control.

- global:

  glob == FALSE -\> only positive correlation; glob == TRUE -\> ±
  partial correlations

## Value

a list object containing two matrics
