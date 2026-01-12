# Create consensus network of resampling-based Markov network

Create consensus network of resampling-based Markov network

## Usage

``` r
consensus_net_ggm(
  ggm_networks,
  CI = 0.95,
  pos_cut = 0,
  neg_cut = 0,
  node_annot = NULL,
  filter = c("pval", "fdr", "none"),
  threshold = 0.05
)
```

## Arguments

- ggm_networks:

  the output from "ensemble_sggm.R"

- CI:

  confidence interval of interests

- pos_cut:

  a threshold for positive partial correlation

- neg_cut:

  a threshold for negative partial correlation

- node_annot:

  a data frame have at least two columns: 'id' and 'symbol', and 'id'
  must match the colnames in "mat"

- filter:

  filter method

- threshold:

  threshold of the selected filter

## Value

a list object including an igraph object as the consensus network
