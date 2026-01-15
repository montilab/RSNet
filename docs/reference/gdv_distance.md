# Find the graphlet degree distance between two vertices

Find the graphlet degree distance between two vertices

## Usage

``` r
gdv_distance(v_gdv, u_gdv, w = 1, similarity = FALSE)
```

## Arguments

- v_gdv:

  the graphlet degree vector of the first node, denoted as v

- u_gdv:

  the graphlet degree vector of the second node, denoted as u

- w:

  the weight vector assigned to each orbit, w = 1 for all orbits by
  default

- similarity:

  if True, return the similarity score, and distance otherwise

## Value

the distance or similarity score

## References

Milenković, Tijana, and Nataša Pržulj. "Uncovering biological network
function via graphlet degree signatures."
