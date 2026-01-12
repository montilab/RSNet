# Simulate network via stochastic block model

Simulate network via stochastic block model

## Usage

``` r
model_sbm(p = 300, m = 6, sizes = NULL, p_in = 0.8, p_out = 0.05)
```

## Arguments

- p:

  The number of vertices in the graph

- m:

  The number of modules

- sizes:

  Size of each module, default = equal-sized modules

- p_in:

  within-module edge density

- p_out:

  between-module edges

## Value

An igraph object
