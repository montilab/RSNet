# Find graphlet degree vector (GDV) of each vertex and graphlet correlation matrix of the given network

Find graphlet degree vector (GDV) of each vertex and graphlet
correlation matrix of the given network

## Usage

``` r
gdvm_gcm(
  ig,
  level = c("4", "5"),
  redundant_orbit = FALSE,
  include_gcm = FALSE,
  correlation = "spearman"
)
```

## Arguments

- ig:

  an igraph object

- level:

  node-level of graphlet

- redundant_orbit:

  if TRUE, keep the redundant orbits, and remove otherwise

- include_gcm:

  if TRUE, compute the GCM

- correlation:

  the correlation measure between each orbit, "spearman" correlation is
  applied by default

## Value

a list of two matrics, gdvm = graphlet degree vector matrix, gcm =
graphlet correlation matrix

## References

Yaveroğlu, Ömer Nebil, et al. "Revealing the hidden language of complex
networks."
