# Compute essential centrality of a network

Compute essential centrality of a network

## Usage

``` r
centrality(ig, weights = NULL)
```

## Arguments

- ig:

  an igraph object

- weights:

  either a numeric edge-weight vector (length = \|E\|) or the name of an
  edge attribute to use as weights (character scalar). Use NULL for
  unweighted.

## Value

a data frame with rows = nodes and columns = centrality scores
