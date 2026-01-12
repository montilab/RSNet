# Resampling-based Markov Network

Resampling-based Markov Network

## Usage

``` r
ensemble_ggm(
  dat,
  num_iteration = 5,
  boot = FALSE,
  sub_ratio = NULL,
  sample_class = NULL,
  correlated = FALSE,
  cluster_ratio = 1,
  estimate_CI = FALSE,
  method = c("D-S_NW_SL", "B_NW_SL"),
  alpha = 0.05,
  global = TRUE,
  n_cores = NULL
)
```

## Arguments

- dat:

  a n x p data frame/matrix with row and column names, n = number of
  samples, and p = number of features

- num_iteration:

  an integer, indicating the number of iteration of resampling

- boot:

  a logical variable, if TRUE, then perform bootstrap resampling, else,
  perform subsampling

- sub_ratio:

  a numerical value between 0 and 1, indicating the subsampling ratio

- sample_class:

  a atomic vector, indicating the class of each sample, if != NULL, then
  stratified sampling is performed

- correlated:

  a boolean variable indicating if samples are correlated

- cluster_ratio:

  the ratio of clusters to be bootstrapped

- estimate_CI:

  if TRUE, save the partial correlation of each edge over the iteration,
  else, return the average partial correlation

- method:

  Methods for statistical inference

- alpha:

  A user-supplied sequence of pre-specified alpha levels for FDR control

- global:

  global == TRUE -\> ± partial correlations; global == FALSE -\> only
  positive correlations

- n_cores:

  number of cores for parallel computing

## Value

a list object
