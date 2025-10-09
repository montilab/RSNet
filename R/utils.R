############################# Structure Inference #############################

#' Perform data resampling
#'
#' @param sample_data a n x 2 data frame is "sample_class" is provided, else, n x 1 data frame of sample IDs
#' @param sample_class a vector indicates the class of each sample, if != NULL, then stratified sampling is performed
#' @param boot a logical variable, if TRUE, then perform bootstrap resampling
#' @param sub_ratio a numerical value between 0 and 1, indicating the subsampling ratio, or NULL


#' @return a vector containing the resampled row.names(sample_data)

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter group_by mutate slice_sample pull

#' @export
resample <- function(sample_data,
                     sample_class,
                     boot=TRUE,
                     sub_ratio=NULL){

  stopifnot(is.data.frame(sample_data))
  if(!is.null(sample_class)){
    stopifnot(is.vector(sample_class))
  }
  if(!is.null(sub_ratio)){
    if( sub_ratio < 0 | sub_ratio > 1) stop(' `sub_ratio` should be between 0 and 1 ')

  }
  if(boot & is.null(sub_ratio)){
    ## bootstrapping, sample n samples w/ replacement
    if(!is.null(sample_class)){
      bagging <- sample_data %>%
        dplyr::group_by(class) %>%
        dplyr::slice_sample(prop=1, replace = TRUE) %>%
        dplyr::pull(var = ID)
    }
    else{
      bagging <- sample_data %>%
        dplyr::slice_sample(prop=1, replace = TRUE) %>%
        dplyr::pull(var = ID)
    }
  }else if(!boot & !is.null(sub_ratio)){
    ## sub-sampling, ratio = sub_ratio, sample (sub_ratio*n) w/o replacement
    if(!is.null(sample_class)){
      bagging <- sample_data %>%
        dplyr::group_by(class) %>%
        dplyr::slice_sample(prop=sub_ratio, replace = FALSE) %>%
        dplyr::pull(var = ID)
    }
    else{
      bagging <- sample_data %>%
        dplyr::slice_sample(prop=sub_ratio, replace = FALSE) %>%
        dplyr::pull(var = ID)
    }
  }else if(boot & !is.null(sub_ratio)){
    ## sample (sub*ratio*n) samples w/ replacement
    if(!is.null(sample_class)){
      bagging <- sample_data %>%
        dplyr::group_by(class) %>%
        dplyr::slice_sample(prop=sub_ratio, replace = TRUE) %>%
        dplyr::pull(var = ID)
    }else{
      bagging <- sample_data %>%
        dplyr::slice_sample(prop=sub_ratio, replace = TRUE) %>%
        dplyr::pull(var = ID)
    }
  }else{
    stop("Invalid sampling method")
  }

  return(bagging)
}

#' Perform data resampling from correlated data
#'
#' @param sample_data a data frame of sample meta-data
#' @param id_col the column name specified the unique sample identifiers
#' @param cluster_col the column name represents the "cluster" of each sample
#' @param boot a logical variable, if TRUE, then perform bootstrap resampling
#' @param ratio the ratio of clusters subject to bootstrapping
#'
#' @return a list object containing bootstrapping cluster and associated samples

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr pull filter rename count inner_join
#' @importFrom tibble as_tibble
#' @importFrom tidyr uncount

#' @export
resample_cluster <- function(sample_data,
                             id_col = "ID",
                             cluster_col = "cluster",
                             boot = TRUE,
                             ratio = NULL){

  stopifnot(is.data.frame(sample_data))
  if(!id_col %in% colnames(sample_data)) stop('Invalid `id_col`')
  if(!cluster_col %in% colnames(sample_data)) stop('Invalid `cluster_col`')

  if(!is.null(ratio)){
    if( ratio < 0 | ratio > 1) stop(' `ratio` should be between 0 and 1 ')
  }

  if(boot & is.null(ratio)){
    ## bootstrapping at cluster-level w/ replacement
    bagging_cluster <- sample_data %>%
      dplyr::pull(!!as.name(cluster_col)) %>%
      sample(size = floor(1*length(unique(.))), replace = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(!!cluster_col := value) %>%
      dplyr::count(!!as.name(cluster_col), name = "counts")
  }else if(!boot & !is.null(ratio)){
    ## select c% of clusters w/o replacement
    bagging_cluster <- sample_data %>%
      dplyr::pull(!!as.name(cluster_col)) %>%
      sample(size = floor(ratio*length(unique(.))), replace = FALSE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(!!cluster_col := value) %>%
      dplyr::count(!!as.name(cluster_col), name = "counts")
  }else if(boot & !is.null(ratio)){
    ## select c% of clusters w/ replacement
    bagging_cluster <- sample_data %>%
      dplyr::pull(!!as.name(cluster_col)) %>%
      sample(size = floor(ratio*length(unique(.))), replace = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(!!cluster_col := value) %>%
      dplyr::count(!!as.name(cluster_col), name = "counts")
  }else{
    stop("Invalid sampling method")
  }

  ## select samples according to the bootstrapped clusters
  bagging <- bagging_cluster %>%
    dplyr::inner_join(sample_data, by = cluster_col) %>%
    tidyr::uncount(weights = counts) %>%
    dplyr::pull(!!as.name(id_col))

  return(list(bagging = bagging,
              bagging_cluster =  bagging_cluster))
}

#' Perform adjustments of p-values on a p x p matrix
#'
#' @param mx_p a p x p matrix, (i,j) represents the p-value of the partial correlation between node i and node j
#'
#' @return a p x p matrix, (i,j) represents the adjusted p-value of the partial correlation between node i and node j

#' @import methods utils

#' @export
matrix_p_adjust <- function( mx_p ) {
  ## initialize
  mx_q <- mx_p
  ## adjust upper triangle
  mx_q[upper.tri(mx_q)] <-
    p.adjust(mx_p[upper.tri(mx_p)], method = "BH")
  ## copy to lower triangle
  mx_q[lower.tri(mx_q)] <-
    t(mx_q)[lower.tri(mx_q)]
  ## some checks
  stopifnot(isTRUE(
    all.equal(
      mx_q[upper.tri(mx_q)],
      t(mx_q)[upper.tri(mx_q)]
    )
  ))
  stopifnot(!isTRUE(all.equal(mx_p,mx_q)))
  return( mx_q)
}

#' Reconstruct the symmetric matrix from upper triangular vector
#'
#' @param upper_tri_values a numeric vector of the upper triangle of the matrix
#' @param variable_names row&column names
#' @param diagl the diagnoal elemenet, a number of a numerci vector

#' @return a symmetric matrix

#' @import methods utils

#' @export
upper_tri_to_matrix <- function(upper_tri_values,
                                variable_names =NULL,
                                diagl=1){
  p <- (1 + sqrt(1 + 8 * length(upper_tri_values))) / 2
  if( (length(diagl)>1) & (length(diagl) != p) ) stop("invalid dignoal!")

  mat <- matrix(0, p, p)

  if(!is.null(variable_names)){
    row.names(mat) <- variable_names
    colnames(mat) <- variable_names
  }

  ## Fill the diagonal
  diag(mat) <- diagl

  ## Fill the upper triangular part
  mat[upper.tri(mat, diag = FALSE)] <- upper_tri_values

  # Fill the lower triangular part (mirror the upper triangular part)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  return(mat)
}


######################## Graphlet ##############################################

#' Pair-wise GDV distance
#'
#' @param gdvm the graphlet degree vector matrix
#' @param w the weight vector
#' @param similarity if True, return the similarity score, and distance otherwise
#'
#' @return a n x n distance/similarity matrix, n = # of node
#' @export
pgdv_distance <- function(gdvm, w=1, similarity=FALSE){

  # Compute pairwise GDV distance/similarity
  mat <- outer(
    1:nrow(gdvm), 1:nrow(gdvm),
    Vectorize(function(i, j) gdv_distance(gdvm[i, ], gdvm[j, ],w=w,similarity = similarity))
  )

  return(mat)

}


######################## Community Detection ###################################

#' Remove isolated nodes from a network
#'
#' @param ig a igraph object
#' @param mode the connection mode
#' @param vertex_symbol the node attribute that represents the symbol of the node
#'
#' @return a list object

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr slice pull
#' @importFrom igraph is.connected degree as_data_frame delete.vertices

#' @export
remove_isolated <- function(ig,
                            mode=c('total','in','out'),
                            vertex_symbol=NULL){

  if(igraph::is.connected(ig)){return(list(ig = ig))}

  if(is.null(vertex_symbol)){
    V(ig)$name <- as.character(seq_len(igraph::vcount(ig)))
    vertex_symbol <- "name"
  }

  if (!is.null(vertex_symbol) && !(vertex_symbol %in% igraph::vertex_attr_names(ig))) {
    stop("'vertex_symbol' NOT in `vertex_attr` of `ig`")
  }

  ## connection mode
  match.arg(mode)

  ## find isolated vertices
  isolated <-  which(igraph::degree(ig, mode=mode)==0)
  ## save the information of isolated vertices
  isolated_IDs <- igraph::as_data_frame(ig, what="vertices") %>%
    dplyr::slice(isolated) %>%
    dplyr::pull(vertex_symbol)
  ## remove isolated vertices
  ig <-  igraph::delete.vertices(ig, isolated)

  return(list(ig = ig, isolated_IDs  = isolated_IDs ))
}


#' Remove small communities in a network
#'
#' @param igs a list containing at least: an igraph object with name "ig", and "ig" must have vertex attribute of community partition
#' @param threshold the minimum community size
#' @param community the name of vertex attribute "community"

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr count filter pull
#' @importFrom igraph as_data_frame delete_vertices

#' @return a list object
#' @export
remove_small_community <- function(igs,
                                   threshold = 5,
                                   community = 'community'){

  large_community <- igraph::as_data_frame(igs$ig, what="vertices") %>%
    dplyr::count(!!as.name(community)) %>%
    dplyr::filter(n > threshold) %>%
    dplyr::pull(!!as.name(community))

  remove_vertices <- igraph::as_data_frame(igs$ig, what="vertices") %>%
    dplyr::filter(!(!!as.name(community) %in% large_community)) %>%
    dplyr::pull(name)

  igs$ig <- igraph::delete_vertices(igs$ig, remove_vertices)

  igs$community_partition <- igs$community_partition[names(igs$community_partition) %in% paste0("C",large_community) == TRUE]


  return(igs)
}


########################### Plot ###############################################

#' An empty ggplot
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 ggplot theme_void
#'
#' @export
ggempty <- function() {
  ggplot() +
    theme_void()
}


######################## Miscellaneous  ########################################

#' Calculate the Jaccard similarity between vectors or matrices
#'
#' @param x a vector or a adjacency matrix
#' @param y a vector or a adjacency matrix

#' @return the Jaccard similarity between x and y

#' @import methods utils

#' @export
jaccard_similarity <- function(x,y){

  if(is.vector(x) & is.vector(y)){
    intersect_norm <- length(intersect(x,y))
    union_norm <- length(union(x,y))
    return(intersect_norm/union_norm)
  }
  if(is.matrix(x) & is.matrix(y)){
    if(!all(dim(x)==dim(y))) stop(" The dimension of matrices must be the same!")
    intersect_size <- sum(x & y)
    union_size <- sum(x | y)

    if (union_size == 0) {
      return(0)  # Handle the case when both matrices are empty
    } else {
      return(intersect_size / union_size)
    }
  }else stop(" 'x' and 'y' must both be 'vector' or 'matrix' !")
}


#' Function to compute Euclidean distance with "pairwise.complete.obs" method
#'
#' @param v1 the first vector
#' @param v2 the second vector
#'
#' @return the Euclidean distance
#' @export
euclidean_pairwise_complete_obs <- function(v1, v2) {
  # Identify valid (non-NA) pairs

  valid_idx <- !is.na(v1) & !is.na(v2)
  # If there are no valid pairs (all are NA), return NA
  if (sum(valid_idx) == 0) {return(NA)}
  # Compute the Euclidean distance on the valid (non-NA) values
  return(sqrt(sum((v1[valid_idx] - v2[valid_idx])^2)))
}






