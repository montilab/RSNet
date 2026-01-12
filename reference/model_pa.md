# Simulate network via preferential attachment model

Simulate network via preferential attachment model

## Usage

``` r
model_pa(p = 300, power = 1, z_appeal = 1, DAG = FALSE)
```

## Arguments

- p:

  The number of vertices in the graph

- power:

  The power of the preferential attachment

- z_appeal:

  The attractiveness of the vertices with no adjacent edges

- DAG:

  If TRUE, simulate an directed acyclic graph, or DAG

## Value

An igraph object
