# Perform data resampling

Perform data resampling

## Usage

``` r
resample(sample_data, sample_class, boot = TRUE, sub_ratio = NULL)
```

## Arguments

- sample_data:

  a n x 2 data frame is "sample_class" is provided, else, n x 1 data
  frame of sample IDs

- sample_class:

  a vector indicates the class of each sample, if != NULL, then
  stratified sampling is performed

- boot:

  a logical variable, if TRUE, then perform bootstrap resampling

- sub_ratio:

  a numerical value between 0 and 1, indicating the subsampling ratio,
  or NULL

## Value

a vector containing the resampled row.names(sample_data)
