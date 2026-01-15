# Create consensus network of resampling-based conditional Gaussian Bayesian network

Create consensus network of resampling-based conditional Gaussian
Bayesian network

## Usage

``` r
consensus_net_cgbn(
  networks,
  reference_network = NULL,
  method = c("all", "average"),
  cut = 0.5
)
```

## Arguments

- networks:

  a list of igraph directed graphs (same vertex set by name)

- reference_network:

  an igraph directed graph (optional, used when method = "all")

- method:

  the method to construct the consensus network: "all" or "average"

- cut:

  minimum frequency (in \[0,1\]) to retain an edge in the consensus

## Value

a list with consensus_network (igraph) and method
