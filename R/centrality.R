#' Compute essential centrality of a network
#'
#' @param ig an igraph object
#' @param weights either a numeric edge-weight vector (length = |E|) or
#'   the name of an edge attribute to use as weights (character scalar). Use NULL for unweighted.
#'
#' @importFrom igraph degree strength eigen_centrality betweenness closeness page_rank
#' @importFrom igraph edge_attr edge_attr_names V
#' @import methods utils
#'
#' @return a data frame with rows = nodes and columns = centrality scores
#' @export
centrality <- function(ig, weights = NULL) {
  if (!inherits(ig, "igraph")) stop("'ig' must be an igraph object")

  ## Resolve weights: allow character attribute name or numeric vector
  w <- weights
  if (is.character(weights)) {
    if (length(weights) != 1L)
      stop("'weights' must be a single attribute name or a numeric vector")
    if (!weights %in% igraph::edge_attr_names(ig))
      stop(sprintf("Edge attribute '%s' not found in graph", weights))
    w <- igraph::edge_attr(ig, weights)
  }

  # degree (unweighted count) and strength (weighted degree)
  degree <- igraph::degree(ig)
  strength <- igraph::strength(ig, weights = w)

  # eigenvector centrality
  eigen <- igraph::eigen_centrality(ig, weights = w)$vector

  # betweenness / closeness: note weights are interpreted as distances (costs)
  betweenness <- igraph::betweenness(ig, weights = w)
  closeness <- igraph::closeness(ig, weights = w)

  # PageRank (non-negative weights expected)
  pagerank <- igraph::page_rank(ig, weights = w)$vector

  df <- data.frame(
    degree      = degree,
    strength    = strength,
    eigen       = eigen,
    betweenness = betweenness,
    closeness   = closeness,
    pagerank    = pagerank,
    check.names = FALSE
  )

  # attach vertex names (if present)
  vn <- igraph::V(ig)$name
  if (!is.null(vn)) rownames(df) <- vn

  return(df)
}
