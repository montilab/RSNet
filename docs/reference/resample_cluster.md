# Perform data resampling from correlated data

Perform data resampling from correlated data

## Usage

``` r
resample_cluster(
  sample_data,
  id_col = "ID",
  cluster_col = "cluster",
  boot = TRUE,
  ratio = NULL
)
```

## Arguments

- sample_data:

  a data frame of sample meta-data

- id_col:

  the column name specified the unique sample identifiers

- cluster_col:

  the column name represents the "cluster" of each sample

- boot:

  a logical variable, if TRUE, then perform bootstrap resampling

- ratio:

  the ratio of clusters subject to bootstrapping

## Value

a list object containing bootstrapping cluster and associated samples
