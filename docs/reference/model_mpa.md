# Simulate network via modular preferential attachment model

Simulate network via modular preferential attachment model

## Usage

``` r
model_mpa(p = 300, m = 6, q_hub = 0.95, m_links = 6, power = 1.7, z_appeal = 1)
```

## Arguments

- p:

  The number of vertices in the graph

- m:

  The number of modules

- q_hub:

  The quantile for degree cutoff for defining modular hubs

- m_links:

  The number of links between modules

- power:

  The power of the preferential attachment

- z_appeal:

  The attractiveness of the vertices with no adjacent edges

## Value

An igraph object
