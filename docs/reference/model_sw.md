# Simulate network via small world model

Simulate network via small world model

## Usage

``` r
model_sw(p = 300, dim = 1, nei = 3, rw_prob = 0.05)
```

## Arguments

- p:

  The number of vertices in the graph

- dim:

  The dimension of the starting lattice

- nei:

  The neighborhood within which the vertices of the lattice will be
  connected

- rw_prob:

  The probability for rewiring

## Value

An igraph object
