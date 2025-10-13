#' Differential connectivity nulls via permutation/bootstrap + GGM
#'
#' Generates a collection of null networks by repeatedly shuffling the data
#' (via permutation or bootstrap) and fitting a sparse Gaussian graphical model
#' (SILGGM) within each label-defined group. Each iteration returns a named list
#' of two igraph objects (one per group). Edge attributes include:
#' \itemize{
#'   \item \code{qval} — BH-adjusted p-values
#'   \item \code{pval} — p-values from precision-matrix inference
#'   \item \code{pcor} — partial correlations
#'   \item \code{sign} — sign of partial correlation (+1 / -1)
#' }
#'
#' @param dat data.frame; p numeric feature columns + one label column.
#' @param inference_method one of \code{c("D-S_NW_SL", "B_NW_SL")}.
#' @param shuffle_method character; one of \code{c("permutation","bootstrap")}.
#' @param shuffle_iter integer; number of null resamples (default 100).
#' @param group_col character; name of the label column (default \code{"phenotype"}).
#' @param balanced logical; only used when \code{shuffle_method = "permutation"}.
#'   If \code{TRUE}, downsample to equal group sizes before permutation.
#' @param filter character; one of \code{c("pval","fdr","none")}.
#' @param threshold numeric; threshold for edge inclusion when \code{filter != "none"}.
#' @param n_cores integer; number of cores parallel computing
#' @param seed integer; random seet
#'
#' @return A list of length \code{shuffle_iter}; each element is a named list
#'   (by group labels) of igraph objects
#'
#' @importFrom parallel mclapply detectCores
#' @export
null_ggm <- function(dat,
                     inference_method = c("D-S_NW_SL", "B_NW_SL"),
                     shuffle_method = c("permutation", "bootstrap"),
                     shuffle_iter = 100,
                     group_col = "phenotype",
                     balanced = FALSE,
                     filter = c("pval", "fdr", "none"),
                     threshold = 0.05,
                     n_cores =  NULL,
                     seed = NULL) {

  stopifnot(is.data.frame(dat))
  stopifnot(group_col %in% names(dat))

  inference_method <- match.arg(inference_method)
  shuffle_method   <- match.arg(shuffle_method)
  filter           <- match.arg(filter)

  cores <- if (is.null(n_cores)) max(1L, parallel::detectCores() - 1L) else n_cores
  if (!is.null(seed)) set.seed(seed)

  run_one <- function(i) {
    this_seed <- if (is.null(seed)) i else seed + i - 1L

    df_list <- if (shuffle_method == "permutation") {
      permutation_diffanal(
        df        = dat,
        group_col = group_col,
        balanced  = balanced,
        seed      = this_seed
      )
    } else {
      bootstrap_diffanal(
        df        = dat,
        group_col = group_col,
        seed      = this_seed
      )
    }

    # Build one igraph per group via internal wrapper
    lapply(
      df_list,
      .diffanal_wrapper,
      method    = inference_method,
      alpha     = 0.05,
      global    = TRUE,
      filter    = filter,
      threshold = threshold
    )
  }

  idx <- seq_len(shuffle_iter)
  return(parallel::mclapply(idx, run_one, mc.cores = cores, mc.set.seed = TRUE))
}

#' Run SILGGM and build an igraph with edge attributes
#'

#' @param dat    n x p numeric matrix/data.frame (rows = samples, cols = features)
#' @param method character; inference method
#' @param alpha  numeric;
#' @param global logical;
#' @param filter one of \code{c("pval","fdr","none")};
#' @param threshold numeric; threshold for edge inclusion.
#'
#' @return igraph object with edge attributes: \code{qval}, \code{pval}, \code{pcor}, \code{sign}
#' @keywords internal
#'
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist V E
#' @importFrom SILGGM SILGGM
.diffanal_wrapper <- function(dat,
                              method,
                              alpha = 0.05,
                              global = TRUE,
                              filter = c("pval", "fdr", "none"),
                              threshold = 0.05) {

  filter <- match.arg(filter)

  # Ensure matrix
  mat <- if (is.matrix(dat)) dat else as.matrix(dat)
  if (is.null(colnames(mat))) colnames(mat) <- paste0("V", seq_len(ncol(mat)))

  # Fit SILGGM
  pcnet <- suppressMessages(
    SILGGM::SILGGM(mat, method = method, alpha = alpha, global = global)
  )

  # Extract matrices
  p_precision <- pcnet$p_precision
  q_precision <- matrix_p_adjust(p_precision)  # provided elsewhere in your package
  pcor        <- pcnet$partialCor

  # Build adjacency by chosen filter (no symmetrization; assumes SILGGM returns symmetric matrices)
  if (filter == "fdr") {
    adj <- (q_precision <= threshold)
    diag(adj) <- FALSE
    A <- adj * 1L
    ig <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  } else if (filter == "pval") {
    adj <- (p_precision <= threshold)
    diag(adj) <- FALSE
    A <- adj * 1L
    ig <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  } else { # filter == "none"
    # Keep weights; by default use p_precision (smaller = stronger). If you prefer
    # partial correlations as weights, replace the next line with: W <- abs(pcor)
    W <- p_precision
    diag(W) <- 0
    ig <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", diag = FALSE, weighted = TRUE)
  }

  igraph::V(ig)$name <- colnames(mat)

  # Attach edge attributes aligned to igraph's edge order
  if (igraph::ecount(ig) > 0L) {
    el <- igraph::as_edgelist(ig, names = FALSE)
    pick <- function(M) M[cbind(el[, 1], el[, 2])]
    igraph::E(ig)$qval <- pick(q_precision)
    igraph::E(ig)$pval <- pick(p_precision)
    igraph::E(ig)$pcor <- pick(pcor)
    igraph::E(ig)$sign <- ifelse(igraph::E(ig)$pcor > 0, 1L, -1L)
  } else {
    igraph::E(ig)$qval <- numeric(0)
    igraph::E(ig)$pval <- numeric(0)
    igraph::E(ig)$pcor <- numeric(0)
    igraph::E(ig)$sign <- integer(0)
  }

  return(ig)
}
