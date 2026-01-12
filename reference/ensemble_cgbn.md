# Resampling-based conditional Gaussian Bayesian Network

Resampling-based conditional Gaussian Bayesian Network

## Usage

``` r
ensemble_cgbn(
  dat,
  discrete_variable = NULL,
  num_iteration = 1,
  boot = FALSE,
  sub_ratio = 0.9,
  sample_class = NULL,
  hugin = FALSE,
  constraints = NULL,
  alpha = 0.01,
  tol = 1e-04,
  maxit = 0,
  n_cores = NULL
)
```

## Arguments

- dat:

  a n x p data frame/matrix with row and column names, n = number of
  samples, and p = number of features

- discrete_variable:

  the columns that represents the discrete_variable, must be a character
  vector

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

- hugin:

  a logical variable, if TRUE, the Hugin object will be included in the
  output

- constraints:

  prior structure knowledge, default = NULL, details -\>
  RHugin::learn.structure

- alpha:

  parameter of PC algorithm for structure learning

- tol:

  parameter of EM algorithm for CPT learning, details -\>
  RHugin::learn.cpt

- maxit:

  parameter of EM algorithm for CPT learning, details -\>
  RHugin::learn.cpt

- n_cores:

  number of cores for parallel computing

## Value

a list object
