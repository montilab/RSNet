# Differential connectivity nulls via permutation/bootstrap + GGM

Generates a collection of null networks by repeatedly shuffling the data
(via permutation or bootstrap) and fitting a sparse Gaussian graphical
model (SILGGM) within each label-defined group. Each iteration returns a
named list of two igraph objects (one per group). Edge attributes
include:

- `qval` — BH-adjusted p-values

- `pval` — p-values from precision-matrix inference

- `pcor` — partial correlations

- `sign` — sign of partial correlation (+1 / -1)

## Usage

``` r
null_ggm(
  dat,
  group_col = "phenotype",
  inference_method = c("D-S_NW_SL", "B_NW_SL"),
  shuffle_method = c("permutation", "bootstrap"),
  shuffle_iter = 100,
  balanced = FALSE,
  filter = c("pval", "fdr", "none"),
  threshold = 0.05,
  n_cores = NULL,
  seed = NULL
)
```

## Arguments

- dat:

  data.frame; p numeric feature columns + one label column.

- group_col:

  character; name of the label column (default `"phenotype"`).

- inference_method:

  one of `c("D-S_NW_SL", "B_NW_SL")`.

- shuffle_method:

  character; one of `c("permutation","bootstrap")`.

- shuffle_iter:

  integer; number of null resamples (default 100).

- balanced:

  logical; only used when `shuffle_method = "permutation"`. If `TRUE`,
  downsample to equal group sizes before permutation.

- filter:

  character; one of `c("pval","fdr","none")`.

- threshold:

  numeric; threshold for edge inclusion when `filter != "none"`.

- n_cores:

  integer; number of cores parallel computing

- seed:

  integer; random seet

## Value

A list of length `shuffle_iter`; each element is a named list (by group
labels) of igraph objects
