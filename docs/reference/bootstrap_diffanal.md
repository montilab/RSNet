# Bootstrap for differential connectivity analysis (pooled bootstrap)

Bootstrap for differential connectivity analysis (pooled bootstrap)

## Usage

``` r
bootstrap_diffanal(df, group_col = "phenotype", max_tries = 100, seed = NULL)
```

## Arguments

- df:

  An input n x (p+1) data frame.

- group_col:

  Column that represents the group/phenotype.

- max_tries:

  Max retries to ensure both classes appear in the bootstrap (default
  100).

- seed:

  Optional random seed.

## Value

A named list of two bootstrapped data frames (one per group), with
`group_col` and the temporary `group_boot` removed.
