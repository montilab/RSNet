# Function to count signed graphlets

Function to count signed graphlets

## Usage

``` r
signed_gdvm_gcm(ig, n_cores = NULL, redundant = FALSE, include_gcm = FALSE)
```

## Arguments

- ig:

  an igraph object

- n_cores:

  number of cores for parrallel computing

- redundant:

  if TRUE, include the redundant orbits

- include_gcm:

  if TRUE, compute GCM

## Value

a list of two matrics, gdvm = graphlet degree vector matrix, gcm =
graphlet correlation matrix

## References

Das, Apratim. Efficient enumeration of small graphlets and orbits.
