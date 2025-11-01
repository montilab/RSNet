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
#' Find the graphlet degree distance between two vertices
#'
#' @param v_gdv the graphlet degree vector of the first node, denoted as v
#' @param u_gdv the graphlet degree vector of the second node, denoted as u
#' @param w the weight vector assigned to each orbit, w = 1 for all orbits by default
#' @param similarity if True, return the similarity score, and distance otherwise

#' @return the distance or similarity score

#' @references Milenković, Tijana, and Nataša Pržulj. "Uncovering biological network function via graphlet degree signatures."

#' @export
gdv_distance <- function(v_gdv, u_gdv, w=1, similarity=FALSE){

  ## check the length
  stopifnot(length(v_gdv) == length(u_gdv))

  ## check the weight vector
  if(length(w) > 1 && (length(w) != length(v_gdv))){stop("Invalid weight vector!")}
  if(length(w) == 1){w <- rep(w, length(v_gdv))}

  ## distance of each orbit
  numerator <- abs(log10(v_gdv+1) - log10(u_gdv+1))
  denominator <- log10(pmax(v_gdv, u_gdv) + 2)
  distances <- numerator/denominator

  if(!similarity){
    return(sum(distances)/sum(w))
  }else{
    return(1 -sum(distances)/sum(w) )
  }

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

#' Differential connectivity bootstrap statistics (bias-corrected)
#'
#' Computes a two-sided p-value and confidence interval for a statistic \code{mi}
#' against its null distribution \code{pi}. Internally uses a bias-corrected
#' mapping on the null via \eqn{z_0} and finds the calibration \eqn{z_A} by
#' minimizing a CI width criterion (Brent's method).
#'
#' @param mi Numeric scalar. The observed measure/statistic.
#' @param pi Numeric vector. The null/sampling distribution (e.g., from
#'   permutation or bootstrap). \code{NA}s are removed.
#' @param CI Numeric in (0,1). Confidence level for the interval. Default 0.95.
#'
#' @return A one-row \code{data.frame} with columns:
#' \itemize{
#'   \item \code{test_dir}: "Up" if \code{mi > 0} else "Down"
#'   \item \code{z}: signed z-score aligned with \code{mi}
#'   \item \code{pval}: two-sided p-value
#'   \item \code{fdr}: Benjamini–Hochberg adjusted p-value (computed on the returned row)
#'   \item \code{ci_low}, \code{ci_high}: confidence interval bounds at level \code{conf}
#' }
#'
#' @details
#' Let \eqn{\hat{F}} be the empirical CDF of \code{pi}. Define
#' \deqn{ z_0 = \mathrm{sign}(-\mathbb{E}[pi - mi]) \cdot \left| \Phi^{-1}\{\hat{F}(mi)\} \right|, }
#' and find \eqn{z_A} that minimizes the CI width induced by the transformed
#' quantiles \eqn{\Phi(2 z_0 \pm z_A)}. The two-sided p-value is
#' \eqn{2 \Phi(-|z_A|)} and the (\code{conf})-level CI uses
#' \eqn{z_\alpha = \Phi^{-1}(\alpha)} with \eqn{\alpha = (1 - \mathrm{conf})/2}.
#'
#'
#' @importFrom stats pnorm qnorm quantile optimize p.adjust
#' @keywords internal
.bootBCStats <- function(mi, pi, CI = 0.95) {
  # --- Validate inputs
  if (!is.numeric(mi) || length(mi) != 1L || is.na(mi)) {
    stop("`mi` must be a single, non-NA numeric value.")
  }
  if (!is.numeric(pi) || length(pi) < 10L) {
    stop("`pi` must be a numeric vector with length >= 10.")
  }
  pi <- pi[is.finite(pi)]
  if (length(pi) < 10L) {
    stop("After removing non-finite values, `pi` has < 10 observations.")
  }
  if (!is.numeric(CI) || length(CI) != 1L || CI <= 0 || CI >= 1) {
    stop("`conf` must be a single numeric in (0, 1).")
  }

  # --- Bias-correction term z0 (clip phat to avoid +/-Inf)
  phat_raw <- mean(pi <= mi)
  eps <- 1 / (length(pi) + 1)                      # continuity correction
  phat <- min(max(phat_raw, eps), 1 - eps)

  ## z0 > 0 -> obs > bootstrap mean
  z0 <- abs(stats::qnorm(phat)) * sign(-mean(pi - mi))

  # --- Find zA by minimizing CI width induced by (2*z0 ± zA)
  obj <- function(zA) .getCI(za = zA, z0 = z0, pi = pi, return_interval = FALSE)
  zA  <- stats::optimize(f = obj, interval = c(-10, 10))$minimum

  # --- Two-sided p-value and signed z aligned with mi
  pval <- 2 * stats::pnorm(-abs(zA))
  z_signed <- abs(stats::qnorm(pval / 2)) * sign(mi)

  # --- Confidence interval at chosen level
  alpha <- (1 - CI) / 2
  z_alpha <- stats::qnorm(alpha)
  ci <- .getCI(za = z_alpha, z0 = z0, pi = pi, return_interval = TRUE)

  res <- data.frame(
    test_dir = ifelse(mi > 0, "Up", "Down"),
    z        = z_signed,
    pval     = pval,
    ci_low   = ci[1],
    ci_high  = ci[2]
  )
  res$fdr <- stats::p.adjust(res$pval, method = "BH")
  return(res)
}

#' @keywords internal
.getCI <- function(za, z0, pi, return_interval = FALSE) {
  q1 <- stats::pnorm(z0 * 2 + za)
  q2 <- stats::pnorm(z0 * 2 - za)
  if (return_interval) {
    stats::quantile(pi, sort(c(q1, q2)), names = FALSE, type = 7)
  } else {
    min(abs(stats::quantile(pi, c(q1, q2), names = FALSE, type = 7)))
  }
}




######################## Miscellaneous  ########################################

#' Silently evaluate an expression while returning its value
#'
#' @description
#' Evaluates an expression, suppressing all console output
#' (print, cat, implicit printing) and all messages, while still
#' returning the normal value of the expression.
#'
#' @param exprs Expression to be evaluated.
#' @param file Optional path to a file where captured output/messages should
#'   be written. If \code{NULL}, they are discarded to the OS null device.
#'
#' @return The normal return value of \code{exprs}.
#' @export
capture_all <- function(exprs, file = NULL) {
  # capture unevaluated user expression
  expr <- substitute(exprs)

  # pick null device
  if (is.null(file)) {
    file <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  }

  # this will hold the *actual* return value
  value <- NULL

  # capture *printed* output
  utils::capture.output(
    # while also suppressing messages
    suppressMessages({
      value <- eval(expr, envir = parent.frame())
    }),
    file = file
  )

  # return the real value (nothing printed)
  return(value)
}


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





