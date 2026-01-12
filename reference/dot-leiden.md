# Perform Leiden community detection algorithm on a network

Perform Leiden community detection algorithm on a network

## Usage

``` r
.leiden(
  igs,
  vertex_symbol = "name",
  resolution = 1,
  objective_function = "modularity",
  n_iterations = 5
)
```

## Arguments

- igs:

  a list containing an igraph object named "ig"

- vertex_symbol:

  the node attribute that represents the symbol of the node

- resolution:

  resolution in Leiden algorithm

- objective_function:

  the objective function in Leiden algorithm

- n_iterations:

  number of iteration in Leiden algorithm

## Value

a list object
