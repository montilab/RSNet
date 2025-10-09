#' Resampling-based Markov Network
#'
#' @param dat a n x p data frame/matrix with row and column names, n = number of samples, and p = number of features
#' @param num_iteration an integer, indicating the number of iteration of resampling
#' @param boot a logical variable, if TRUE, then perform bootstrap resampling, else, perform subsampling
#' @param sub_ratio a numerical value between 0 and 1, indicating the subsampling ratio
#' @param sample_class a atomic vector, indicating the class of each sample, if != NULL, then stratified sampling is performed
#' @param correlated a boolean variable indicating if samples are correlated
#' @param cluster_ratio the ratio of clusters to be bootstrapped
#' @param estimate_CI if TRUE, save the partial correlation of each edge over the iteration, else, return the average partial correlation
#' @param method Methods for statistical inference
#' @param alpha A user-supplied sequence of pre-specified alpha levels for FDR control
#' @param global global == TRUE -> ± partial correlations; global == FALSE -> only positive correlations
#' @param n_cores number of cores for parallel computing
#'
#' @return a list object
#'
#' @importFrom magrittr %>%
#' @importFrom SILGGM SILGGM
#' @importFrom parallel mclapply detectCores
#' @export
ensemble_ggm <- function(dat,
                         num_iteration = 5,
                         boot = FALSE,
                         sub_ratio = NULL,
                         sample_class = NULL,
                         correlated = FALSE,
                         cluster_ratio = 1,
                         estimate_CI = FALSE,
                         method = c("D-S_NW_SL", "B_NW_SL"),
                         alpha = 0.05,
                         global = TRUE,
                         n_cores = NULL) {

  ## Input checks
  if (!(is.data.frame(dat) || is.matrix(dat))) stop("'data' must either be a data.frame or matrix")
  if (is.null(colnames(dat))) stop("'data' must have column names")
  if (is.null(rownames(dat))) stop("'data' must have row names")
  if (!is.null(sample_class) && !is.atomic(sample_class)) stop("'sample_class' must be an atomic vector")

  ## Validate method
  method <- match.arg(method)

  ## Prepare parameters
  n <- nrow(dat)
  p <- ncol(dat)
  ecount <- (p * (p - 1)) / 2
  vids <- colnames(dat)

  ## Prepare input sample_data
  if (!is.null(sample_class) && !correlated) {
    sample_data <- data.frame(ID = rownames(dat), class = as.factor(sample_class))
  } else if (!is.null(sample_class) && correlated) {
    sample_data <- data.frame(ID = rownames(dat), cluster = sample_class)
  } else {
    sample_data <- data.frame(ID = rownames(dat))
  }

  ## Define a single iteration
  run_iteration <- function(i) {
    bagging <- NULL
    bagging_cluster <- NULL

    if (correlated) {
      output <- resample_cluster(sample_data = sample_data,
                                 id_col = "ID",
                                 cluster_col = "cluster",
                                 boot = boot,
                                 ratio = cluster_ratio)
      bagging <- output$bagging
      bagging_cluster <- output$bagging_cluster
    } else {
      bagging <- resample(sample_data = sample_data,
                          sample_class = sample_class,
                          boot = boot,
                          sub_ratio = sub_ratio)
    }

    resampled_data <- dat[bagging, , drop = FALSE]

    pcnet <- .sggm_wrapper(mat = resampled_data,
                           method = method,
                           alpha = alpha,
                           global = global)

    return(list(
      bagging = bagging,
      bagging_cluster = bagging_cluster,
      partialCor = pcnet$partialCor,
      z_score_partialCor = pcnet$z_score_partialCor
    ))
  }

  ## Initialize number of cores
  cores <- if (is.null(n_cores)) max(1, parallel::detectCores() - 1) else n_cores

  ## Run iterations in parallel
  results <- parallel::mclapply(1:num_iteration, run_iteration, mc.cores = cores)

  ## Organize results
  partialCor_mat <- matrix(0, nrow = num_iteration, ncol = ecount)
  z_score_partialCor <- numeric(ecount)
  baggings <- list()

  if (correlated){
    bagging_clusters <- list()
  }

  for (i in seq_along(results)) {

    partialCor_mat[i, ] <- results[[i]]$partialCor

    z_score_partialCor <- z_score_partialCor + results[[i]]$z_score_partialCor

    baggings[[paste0("iter_", i)]] <- results[[i]]$bagging

    if (correlated) {
      bagging_clusters[[paste0("iter_", i)]] <- results[[i]]$bagging_cluster
    }
  }


  ## Construct output
  if (estimate_CI) {
    res <- list(
      method = method,
      estimate_CI = TRUE,
      n = n,
      p = p,
      partialCor_mat = partialCor_mat,
      avg_z_score_partialCor = z_score_partialCor / num_iteration,
      vids = vids,
      baggings = baggings
    )
  } else {
    res <- list(
      method = method,
      estimate_CI = FALSE,
      n = n,
      p = p,
      avg_partialCor = colMeans(partialCor_mat),
      avg_z_score_partialCor = z_score_partialCor / num_iteration,
      vids = vids,
      baggings = baggings
    )
  }

  if (correlated) {
    res$bagging_clusters <- bagging_clusters
  }

  return(res)
}


#' A Wrapper function of SILGGM
#'
#' @param mat a n x p matrix, n = number of samples, and p = number of features
#' @param method the inference method
#' @param alpha A user-supplied sequence of pre-sepecified alpha levels for FDR control.
#' @param global glob == FALSE -> only positive correlation; glob == TRUE -> ± partial correlations
#'
#' @return a list object containing two matrics
#' @keywords internal
.sggm_wrapper <- function(mat,
                          method,
                          alpha = 0.05,
                          global = TRUE) {

  if(!is.matrix(mat)){mat <- as.matrix(mat)}
  n <- nrow(mat)
  output <- suppressMessages(SILGGM::SILGGM(mat, method = method, alpha = alpha, global = global))

  if(method == "D-S_NW_SL"){
    partialCor <- .upper_tri_vec(output$partialCor)
    z_score_partialCor <- sapply(partialCor, .pcor_zscore, n)
    pcnet <- list(partialCor = partialCor,
                  z_score_partialCor = z_score_partialCor)
  }
  if(method == "B_NW_SL"){
    partialCor <- .upper_tri_vec(output$partialCor)
    z_score_partialCor <- .upper_tri_vec(output$z_score_partialCor)
    pcnet <- list(partialCor = partialCor,
                  z_score_partialCor = z_score_partialCor)
  }

  return(pcnet)
}

#' Raise of number to its n-th power
#'
#' @param x the input number
#' @param n the order of power
#'
#' @return a number
#' @keywords internal
.pow <- function(x,n){
  return(x^n)
}

#' Compute the z-score of a partial correlation, see https://github.com/cran/SILGGM/blob/master/src/SILGGMCpp.cpp
#'
#' @param pcor the partial correlation
#' @param n the number of samples
#'
#' @return the z-score of the corresponding partial correlation
#' @keywords internal
.pcor_zscore <- function(pcor,n){
  std_new <- sqrt(.pow((1-.pow(pcor,2)),2)/n)
  return(pcor/std_new)
}

#' Take the upper triangular vector of a symmetric matrix
#'
#' @param mat the input symmetric matrix
#'
#' @return a numeric vector
#' @keywords internal
.upper_tri_vec <- function(mat){
  vec <- mat[upper.tri(mat)]
  if(length(vec) != (ncol(mat)*(ncol(mat)-1))/2) stop ("invalid length!")
  return(vec)
}
