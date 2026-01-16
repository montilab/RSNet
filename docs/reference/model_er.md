# Simulate network via erdos-renyi model

Simulate network via erdos-renyi model

## Usage

``` r
model_er(p = 300, e_prob = 0.01, DAG = FALSE)
```

## Arguments

- p:

  The number of vertices in the graph

- e_prob:

  A probability for drawing an edge between two arbitrary vertices

- DAG:

  If TRUE, simulate an directed acyclic graph, or DAG

## Value

An igraph object
