# Plot Consensus Network

Visualizes a consensus network using visNetwork, with optional
query-based subsetting and optional confidence-interval display for
partial correlations.

## Usage

``` r
plot_cn(
  ig,
  query = NULL,
  order = 1,
  vertex_symbol = "name",
  edge_label = "pcor",
  edge_width = NULL,
  CI_show = FALSE,
  main = "Title",
  query_color = "#FF9999",
  neighbor_color = "#56B4E9",
  all_color = "#56B4E9",
  hover = TRUE,
  stabilization = FALSE,
  igraph_layout = "layout_with_fr",
  layout_seed = 42
)
```

## Arguments

- ig:

  an igraph object

- query:

  character vector of node names/symbols to center the subgraph on

- order:

  neighborhood order around the query nodes

- vertex_symbol:

  vertex attribute to use as node label/symbol

- edge_label:

  edge attribute to show as label; default is "pcor"

- edge_width:

  optional edge attribute for edge thickness

- CI_show:

  logical; if TRUE, show CI \*\*only\*\* when edge_label = "pcor" and
  both "lower_quantile" and "upper_quantile" are available

- main:

  title

- query_color:

  color for queried nodes

- neighbor_color:

  color for neighbor nodes

- all_color:

  color for all nodes when no query is given

- hover:

  logical

- stabilization:

  logical

- igraph_layout:

  layout name to pass to visIgraphLayout

- layout_seed:

  random seed for layout

## Value

list with p (visNetwork object) and data (nodes/edges)
