#' Distance between two Graphlet Correlation Matrices
#'
#' @param gcm1 the GCM of a network
#' @param gcm2 the GCM of another network
#'
#' @returns numeric; the Euclidean distance between the upper triangular of gcm1 and gcm2
#' @export
gcm_distance <- function(gcm1, gcm2){

  stopifnot(is.matrix(gcm1), is.matrix(gcm2))
  stopifnot(is.numeric(gcm1), is.numeric(gcm2))
  stopifnot(identical(dim(gcm1), dim(gcm2)))
  stopifnot(nrow(gcm1) == ncol(gcm1))

  v1 <- gcm1[upper.tri(gcm1)]
  v2 <- gcm1[upper.tri(gcm2)]

  return(.euclidean_pairwise_complete_obs(v1, v2))
}



#' Function to compute Euclidean distance with "pairwise.complete.obs" method
#'
#' @param v1 the first vector
#' @param v2 the second vector
#'
#' @return the Euclidean distance
#' @export
.euclidean_pairwise_complete_obs <- function(v1, v2) {
  # Identify valid (non-NA) pairs

  valid_idx <- !is.na(v1) & !is.na(v2)
  # If there are no valid pairs (all are NA), return NA
  if (sum(valid_idx) == 0) {return(NA)}
  # Compute the Euclidean distance on the valid (non-NA) values
  return(sqrt(sum((v1[valid_idx] - v2[valid_idx])^2)))
}
