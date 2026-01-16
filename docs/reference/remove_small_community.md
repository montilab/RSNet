# Remove small communities in a network

Remove small communities in a network

## Usage

``` r
remove_small_community(igs, threshold = 5, community = "community")
```

## Arguments

- igs:

  a list containing at least: an igraph object with name "ig", and "ig"
  must have vertex attribute of community partition

- threshold:

  the minimum community size

- community:

  the name of vertex attribute "community"

## Value

a list object
