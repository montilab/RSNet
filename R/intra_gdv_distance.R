#' Pair-wise GDV distance
#'
#' @param gdvm the graphlet degree vector matrix
#' @param w the weight vector
#' @param similarity if True, return the similarity score, and distance otherwise
#'
#' @return a n x n distance/similarity matrix, n = # of node
#' @export
intra_gdv_distance <- function(gdvm, w=1, similarity=FALSE){

  # Compute pairwise GDV distance/similarity
  mat <- outer(
    1:nrow(gdvm), 1:nrow(gdvm),
    Vectorize(function(i, j) gdv_distance(gdvm[i, ], gdvm[j, ],w=w,similarity = similarity))
  )

  return(mat)

}
