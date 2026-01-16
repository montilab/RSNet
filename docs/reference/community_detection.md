# Perform community detection on a network

Perform community detection on a network

## Usage

``` r
community_detection(
  ig,
  mode = c("total", "in", "out"),
  vertex_symbol = "name",
  method = c("walk_trap", "louvain", "leiden"),
  small_community = NULL,
  steps = 10,
  resolution = 1,
  n_iterations = 5
)
```

## Arguments

- ig:

  an igraph object

- mode:

  the connection mode used to defined "isolated" nodes

- vertex_symbol:

  the node attribute that represents the symbol of the node

- method:

  the community detection method

- small_community:

  the threshold of small communities to be removed, is NULL, small
  communities will not be removed

- steps:

  the number of steps in walk-trap algorithm

- resolution:

  the resolution in Louvain and/or Leiden algorithm

- n_iterations:

  the number of iteration in Leiden algorithm

## Value

a list object
