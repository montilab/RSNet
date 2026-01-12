# Pair-wise GDV distance

Pair-wise GDV distance

## Usage

``` r
intra_gdv_distance(gdvm, w = 1, similarity = FALSE)
```

## Arguments

- gdvm:

  the graphlet degree vector matrix

- w:

  the weight vector

- similarity:

  if True, return the similarity score, and distance otherwise

## Value

a n x n distance/similarity matrix, n = \# of node
