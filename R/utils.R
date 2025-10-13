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

#' Distinct colors for large categorical datasets
#'
#' Colors derived from the R package colorspace
#' colorspace: A Toolbox for Manipulating and Assessing Colors and Palettes
#' https://cran.r-project.org/web/packages/colorspace/index.html
#' License: BSD_3_clause + https://cran.r-project.org/web/packages/colorspace/LICENSE
#'
#' @param n Number of colors desired
#'
#' @return Color palette
#'
#' @export
distinct_colors <- function(n) {
  x <- c("#353E7C","#007094","#009B95","#00BE7D","#96D84B","#FDE333","#040404","#3E134F","#851170","#C53270",
         "#F36E35","#F8B83C","#A23D47","#A96C00","#9A9800","#75C165","#50E2BB","#B0F4FA","#26185F","#005D9E",
         "#18BDB0","#9ADCBB","#D7F4CF","#DD008F","#DB6AC0","#D5A6DB","#F6F6FC","#CBA079","#928261","#605F4C",
         "#363932","#001889","#87008D","#DAFF47","#88002D","#FE945C","#FFE2C0","#004533","#006F69","#0091AD",
         "#EDD788","#AB4A3D","#73243C","#AC0535","#EB2C31","#EF4868","#F56553","#404E9C","#3D7CB8","#4BA5BF",
         "#55C7B1","#A3E292","#FAEF8B","#0A1230","#3C2871","#7A3392","#B7509A","#E68375","#F1C687","#985277",
         "#A37B49","#9BA453","#86CBA0","#7DEBEA","#C0FCFC","#2C2C7D","#396ABC","#5DC6DB","#AAE6EA","#DEFDFD",
         "#CA40B4","#CF82E3","#D1B6F3","#F7FFFF","#C7AEAE","#928F93","#626C7B","#3A465F","#25309C","#7B31A9",
         "#DEFD99","#7D245B","#F3A698","#FDF0F3","#215061","#367A94","#4B9CD2","#EBE4C4","#A15D70","#6B3868")
  return(rep(x, floor(n/length(x))+1)[1:n])
}

#' Normalize values between a given range
#'
#' @param x Values to normalize
#' @param a Min of range
#' @param b max of range
#'
#' @return Normalized values
#'
#' @export
normalize_range <- function(x, a=0, b=1) {
  (b-a)*( (x-min(x)) / (max(x)-min(x)) )+a
}


#' Colorize numerical values
#'
#' @param values values to be colorized
#' @param resolution Limit resolution for small values
#'
#' @return Colorized values
#' @importFrom grDevices heat.colors
#' @export
colorize <- function(values, resolution=4) {
  multiplier <- 100*resolution
  colors <- rev(heat.colors(multiplier+1))
  colors[round(normalize_range(values)*multiplier, 0)+1]
}


#' Assign unique color to each community
#'
#' @param ig an igraph object
#' @param community_attr_name the node attribute that represents the community partition

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr select pull left_join rename
#' @importFrom igraph vertex_attr_names as_data_frame V

#' @return an igraph object with "color" node attribute
#' @export
colorize_community <- function(ig,
                               community_attr_name = "community"){

  stopifnot(is(ig,"igraph"))
  if(!(community_attr_name %in% igraph::vertex_attr_names(ig))) stop("`community_attr_name` should be in `vertex_attr` of `ig`")

  ## obtain unique community labels
  communities_unique <- igraph::as_data_frame(ig, what='vertices') %>%
    dplyr::pull(!!as.name(community_attr_name)) %>%
    unique()

  ## create color palette
  color_pal <- data.frame(community = communities_unique,
                          color = distinct_colors(length(communities_unique)))

  igraph::V(ig)$color <- igraph::as_data_frame(ig, what='vertices') %>%
    dplyr::rename(community = !!as.name(community_attr_name)) %>%
    dplyr::select(community) %>%
    dplyr::left_join(color_pal, by="community") %>%
    dplyr::pull(color)

  return(ig)

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

#################### Differential Connectivity  ###############################

#' Permutation for differential connectivity analysis
#'
#' @param df an input n by (p+1) data frame
#' @param group_col the column that represents the "group/phenotype" to be shuffled
#' @param balanced a logical parameter for regular or balanced permutation
#' @param seed random seed
#'
#' @returns an list object containing two shuffled data frame for each group
#' @export
permutation_diffanal <- function(df,
                                  group_col = "phenotype",
                                  balanced = FALSE,
                                  seed = NULL) {

  stopifnot(is.data.frame(df))
  stopifnot(group_col %in% names(df))
  if (!is.null(seed)) set.seed(seed)

  # Check: all columns except group_col are numeric
  numeric_check <- vapply(df[ , !(names(df) %in% group_col), drop = FALSE],
                          is.numeric, logical(1))
  if (!all(numeric_check)) {
    bad_cols <- names(numeric_check)[!numeric_check]
    stop(sprintf("All columns except '%s' must be numeric. Non-numeric column(s) found: %s",
                 group_col, paste(bad_cols, collapse = ", ")))
  }

  group <- df[[group_col]]

  if (!is.character(group) && !is.factor(group) &&
      !is.integer(group) && !is.logical(group)) {
    stop("'group_col' must be character, factor, integer, or logical")
  }

  group_vals <- unique(group)
  if (length(group_vals) != 2L) {
    stop(sprintf("`%s` must have exactly 2 unique values (found %d).",
                 group_col, length(group_vals)))
  }

  n1 <- sum(group == group_vals[1])
  n2 <- sum(group == group_vals[2])

  if (!balanced) {
    # Unbalanced: standard label shuffle (preserves counts)
    df$group_perm <- sample(group, length(group), replace = FALSE)
  } else {
    # Balanced: down-sample to equal group sizes
    n_min <- min(n1, n2)
    idx1 <- sample(which(group == group_vals[1]), n_min)
    idx2 <- sample(which(group == group_vals[2]), n_min)
    idx_bal <- c(idx1, idx2)

    df <- df[idx_bal, , drop = FALSE]
    df$group_perm <- sample(df[[group_col]], length(df[[group_col]]), replace = FALSE)
  }

  # Build list by permuted label
  df_list <- lapply(group_vals, function(v) {
    df_sub <- df[df$group_perm == v, , drop = FALSE]
    # Drop both group_col and group_perm columns
    df_sub <- df_sub[, !(names(df_sub) %in% c(group_col, "group_perm")), drop = FALSE]
    return(df_sub)
  })

  names(df_list) <- group_vals
  return(df_list)
}

#' Bootstrap for differential connectivity analysis (pooled bootstrap)
#'
#' @param df        An input n x (p+1) data frame.
#' @param group_col Column that represents the group/phenotype.
#' @param max_tries Max retries to ensure both classes appear in the bootstrap (default 100).
#' @param seed      Optional random seed.
#'
#' @return A named list of two bootstrapped data frames (one per group),
#'         with \code{group_col} and the temporary \code{group_boot} removed.
#' @export
bootstrap_diffanal <- function(df,
                               group_col = "phenotype",
                               max_tries = 100,
                               seed = NULL) {
  stopifnot(is.data.frame(df))
  stopifnot(group_col %in% names(df))
  if (!is.null(seed)) set.seed(seed)

  # 1) Check: all columns except group_col are numeric
  numeric_check <- vapply(df[, !(names(df) %in% group_col), drop = FALSE],
                          is.numeric, logical(1))
  if (!all(numeric_check)) {
    bad_cols <- names(numeric_check)[!numeric_check]
    stop(sprintf("All columns except '%s' must be numeric. Non-numeric column(s) found: %s",
                 group_col, paste(bad_cols, collapse = ", ")))
  }

  # 2) Validate group column
  group <- df[[group_col]]
  if (!is.character(group) && !is.factor(group) &&
      !is.integer(group) && !is.logical(group)) {
    stop("'group_col' must be character, factor, integer, or logical")
  }
  if (anyNA(group)) {
    stop(sprintf("`%s` contains NA; please handle missing labels first.", group_col))
  }

  # Stable naming
  group_vals <- unique(as.character(group))
  if (length(group_vals) != 2L) {
    stop(sprintf("`%s` must have exactly 2 unique values (found %d).",
                 group_col, length(group_vals)))
  }

  # 3) Pooled bootstrap of rows (labels come along; counts vary)
  n <- nrow(df)
  tries <- 0L
  repeat {
    idx_boot <- sample.int(n, size = n, replace = TRUE)
    df_boot  <- df[idx_boot, , drop = FALSE]

    # Inherit labels from resampled rows (no shuffling)
    df_boot$group_boot <- as.character(df_boot[[group_col]])

    tab <- table(df_boot$group_boot)
    # proceed if both target groups are present and non-zero
    if (all(group_vals %in% names(tab)) && all(tab[group_vals] > 0)) break

    tries <- tries + 1L
    if (tries >= max_tries) {
      stop("Bootstrap repeatedly produced a single-class sample; increase 'max_tries' or sample size.")
    }
  }

  # 4) Build list by bootstrapped label; drop label columns
  df_list <- lapply(group_vals, function(v) {
    df_sub <- df_boot[df_boot$group_boot == v, , drop = FALSE]
    df_sub <- df_sub[, !(names(df_sub) %in% c(group_col, "group_boot")), drop = FALSE]
    df_sub
  })

  names(df_list) <- group_vals

  return(df_list)
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






