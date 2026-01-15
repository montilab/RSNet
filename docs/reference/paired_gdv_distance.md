# Compute Paired Graphlet Degree Vector (GDV) Distances

Computes node-wise GDV distances between two Graphlet Degree Vector
Matrices (GDVMs) corresponding to the same set of nodes. The function
aligns rows by node names and returns a named vector of GDV distances.

## Usage

``` r
paired_gdv_distance(gdvm1, gdvm2)
```

## Arguments

- gdvm1:

  A matrix or data frame representing the first GDVM. Rows correspond to
  nodes, and columns to graphlet orbit counts or features.

- gdvm2:

  A matrix or data frame representing the second GDVM. Must have at
  least one overlapping row name with `gdvm1`.

## Value

A named numeric vector giving the GDV distance for each node present in
both GDVMs.

## Details

The function matches nodes by their row names and computes the GDV
distance using
[`gdv_distance`](https://montilab.github.io/RSNet/reference/gdv_distance.md)
for each shared node. Nodes that do not appear in both matrices are
automatically excluded.

## Examples

``` r
if (FALSE) { # \dontrun{
gdv1 <- matrix(runif(50), nrow = 10)
gdv2 <- matrix(runif(50), nrow = 10)
rownames(gdv1) <- rownames(gdv2) <- paste0("node", 1:10)
paired_gdv_dist(gdvm1 = gdv1, gdvm2 = gdv2)
} # }
```
