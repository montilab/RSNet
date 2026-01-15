# Permutation for differential connectivity analysis

Permutation for differential connectivity analysis

## Usage

``` r
permutation_diffanal(
  df,
  group_col = "phenotype",
  balanced = FALSE,
  seed = NULL
)
```

## Arguments

- df:

  an input n by (p+1) data frame

- group_col:

  the column that represents the "group/phenotype" to be shuffled

- balanced:

  a logical parameter for regular or balanced permutation

- seed:

  random seed

## Value

an list object containing two shuffled data frame for each group
