# Run SILGGM and build an igraph with edge attributes

Run SILGGM and build an igraph with edge attributes

## Usage

``` r
.diffanal_wrapper(
  dat,
  method,
  alpha = 0.05,
  global = TRUE,
  filter = c("pval", "fdr", "none"),
  threshold = 0.05
)
```

## Arguments

- dat:

  n x p numeric matrix/data.frame (rows = samples, cols = features)

- method:

  character; inference method

- alpha:

  numeric;

- global:

  logical;

- filter:

  one of `c("pval","fdr","none")`;

- threshold:

  numeric; threshold for edge inclusion.

## Value

igraph object with edge attributes: `qval`, `pval`, `pcor`, `sign`
