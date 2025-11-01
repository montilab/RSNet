#' Compute Paired Graphlet Degree Vector (GDV) Distances
#'
#' Computes node-wise GDV distances between two Graphlet Degree Vector Matrices
#' (GDVMs) corresponding to the same set of nodes. The function aligns rows by
#' node names and returns a named vector of GDV distances.
#'
#' @param gdvm1 A matrix or data frame representing the first GDVM.
#'   Rows correspond to nodes, and columns to graphlet orbit counts or features.
#' @param gdvm2 A matrix or data frame representing the second GDVM.
#'   Must have at least one overlapping row name with \code{gdvm1}.
#' @return A named numeric vector giving the GDV distance for each node present
#'   in both GDVMs.
#' @details
#' The function matches nodes by their row names and computes the GDV distance
#' using \code{\link{gdv_distance}} for each shared node. Nodes that do not
#' appear in both matrices are automatically excluded.
#'
#' @examples
#' \dontrun{
#' gdv1 <- matrix(runif(50), nrow = 10)
#' gdv2 <- matrix(runif(50), nrow = 10)
#' rownames(gdv1) <- rownames(gdv2) <- paste0("node", 1:10)
#' paired_gdv_dist(gdvm1 = gdv1, gdvm2 = gdv2)
#' }
#' @export
paired_gdv_distance <- function(gdvm1, gdvm2) {

  # --- Input validation -------------------------------------------------------

  if (!is.matrix(gdvm1) || !is.matrix(gdvm2)) {
    stop("'gdvm1' and 'gdvm2' must both be matrices or data frames.", call. = FALSE)
  }

  if (is.null(rownames(gdvm1)) || is.null(rownames(gdvm2))) {
    stop("Both 'gdvm1' and 'gdvm2' must have row names representing node IDs.", call. = FALSE)
  }

  # --- Align and match nodes --------------------------------------------------
  inter_rows <- intersect(rownames(gdvm1), rownames(gdvm2))
  if (length(inter_rows) == 0L) {
    stop("No overlapping row names found between the two GDVMs.", call. = FALSE)
  }

  gdvm1 <- gdvm1[inter_rows, , drop = FALSE]
  gdvm2 <- gdvm2[inter_rows, , drop = FALSE]

  # --- Compute GDV distances --------------------------------------------------
  dist_vec <- vapply(
    inter_rows,
    function(node) {
      v1 <- gdvm1[node, , drop = TRUE]
      v2 <- gdvm2[node, , drop = TRUE]
      gdv_distance(v1, v2, w = 1, similarity = FALSE)
    },
    numeric(1)
  )

  # --- Return named vector ----------------------------------------------------
  names(dist_vec) <- inter_rows
  return(dist_vec)
}
