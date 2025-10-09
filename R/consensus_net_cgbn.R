#' Create consensus network of resampling-based conditional Gaussian Bayesian network
#'
#' @param networks a list of igraph directed graphs (same vertex set by name)
#' @param reference_network an igraph directed graph (optional, used when method = "all")
#' @param method the method to construct the consensus network: "all" or "average"
#' @param cut minimum frequency (in [0,1]) to retain an edge in the consensus
#'
#' @import methods utils
#' @importFrom igraph V as_adjacency_matrix graph_from_adjacency_matrix intersection difference
#' @return a list with consensus_network (igraph) and method
#' @export
consensus_net_cgbn <- function(networks,
                               reference_network = NULL,
                               method = c("all", "average"),
                               cut = 0.5) {

  method <- match.arg(method)

  if (!is.list(networks) || length(networks) == 0)
    stop("'networks' must be a non-empty list of igraph objects")

  if (!inherits(networks[[1]], "igraph"))
    stop("'networks[[1]]' must be an igraph object")

  # master vertex order (by name) from first network
  labels <- igraph::V(networks[[1]])$name
  if (is.null(labels)) stop("Vertex names are required in 'networks[[1]]'")

  # sanity: all networks must contain same vertex set (by name)
  for (g in networks) {
    if (!inherits(g, "igraph")) stop("All elements of 'networks' must be igraph graphs")
    ln <- igraph::V(g)$name
    if (is.null(ln) || !setequal(ln, labels))
      stop("All networks must have the same vertex names")
  }

  if (!is.null(reference_network) && !inherits(reference_network, "igraph")) {
    # allow passing a matrix as reference; coerce if needed
    if (is.matrix(reference_network)) {
      reference_network <- igraph::graph_from_adjacency_matrix(reference_network,
                                                               mode = "directed",
                                                               weighted = NULL)
    } else {
      stop("'reference_network' must be an igraph (or adjacency matrix).")
    }
  }

  if (method == "all") {
    consensus_network <- .consensus_all(reference_network, networks, labels, cut)
  } else {
    consensus_network <- .consensus_average(networks, labels, cut)
  }

  return (list(consensus_network = consensus_network, method = method))
}

# ---- helpers -------------------------------------------------------------

# Intersection-based consensus of directed edges against a reference
.consensus_all <- function(reference_network, networks, labels, cut) {
  if (is.null(reference_network))
    stop("When method = 'all', 'reference_network' must be provided")

  # Ensure reference has same vertex set & order
  ref_names <- igraph::V(reference_network)$name
  if (is.null(ref_names) || !setequal(ref_names, labels))
    stop("'reference_network' must have the same vertex names as 'networks'")

  # Keep only edges present in reference; average their presence across networks
  mats <- lapply(networks, function(g) {
    # keep only reference edges
    gi <- igraph::intersection(g, reference_network, byname = TRUE, keep.all.vertices = TRUE)
    A  <- as.matrix(igraph::as_adjacency_matrix(gi, sparse = FALSE))
    # reorder to master labels
    A[labels, labels, drop = FALSE]
  })

  avgA <- Reduce(`+`, mats) / length(mats)
  avgA[avgA < cut] <- 0

  return(igraph::graph_from_adjacency_matrix(avgA, mode = "directed", weighted = TRUE))
}

# Markov-blanket-based consensus (undirected) without MXM/bnlearn
.consensus_average <- function(networks, labels, cut) {
  p <- length(labels)
  CMB <- matrix(0, nrow = p, ncol = p, dimnames = list(labels, labels))

  # For each network, compute MB for every node and accumulate
  for (g in networks) {
    A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    # reorder to master labels
    A <- A[labels, labels, drop = FALSE]
    # binarize
    A[A != 0] <- 1

    for (i in seq_len(p)) {
      mb_idx <- .mb_from_adj(A, i)
      if (length(mb_idx)) {
        CMB[i, mb_idx] <- CMB[i, mb_idx] + 1
        # (optionally symmetrize as you go)
        # CMB[mb_idx, i] <- CMB[mb_idx, i] + 1
      }
    }
  }

  # Average across networks, threshold, return undirected weighted consensus
  CMB <- CMB / length(networks)
  CMB[CMB < cut] <- 0

  return(igraph::graph_from_adjacency_matrix(CMB, mode = "undirected", weighted = TRUE))
}

# Markov blanket of node 'i' from a 0/1 directed adjacency matrix A (rows: from, cols: to)
# Parents: A[, i] == 1
# Children: A[i, ] == 1
# Spouses: co-parents of children, excluding i/parents/children
.mb_from_adj <- function(A, i) {
  parents  <- which(A[, i] != 0)
  children <- which(A[i, ] != 0)

  spouses <- integer(0)
  if (length(children)) {
    for (c in children) {
      spouses <- c(spouses, which(A[, c] != 0))
    }
  }
  spouses <- setdiff(unique(spouses), c(i, parents, children))

  return(sort(unique(c(parents, children, spouses))))
}
