#' Resampling-based conditional Gaussian Bayesian Network
#'
#' @param dat a n x p data frame/matrix with row and column names, n = number of samples, and p = number of features
#' @param discrete_variable the columns that represents the discrete_variable, must be a character vector
#' @param num_iteration an integer, indicating the number of iteration of resampling
#' @param boot a logical variable, if TRUE, then perform bootstrap resampling, else, perform subsampling
#' @param sub_ratio a numerical value between 0 and 1, indicating the subsampling ratio
#' @param sample_class a atomic vector, indicating the class of each sample, if != NULL, then stratified sampling is performed
#' @param hugin a logical variable, if TRUE, the Hugin object will be included in the output
#' @param constraints prior structure knowledge, default = NULL, details -> RHugin::learn.structure
#' @param alpha parameter of PC algorithm for structure learning
#' @param tol parameter of EM algorithm for CPT learning, details -> RHugin::learn.cpt
#' @param maxit parameter of EM algorithm for CPT learning, details -> RHugin::learn.cpt
#' @param n_cores number of cores for parallel computing
#'
#' @return a list object
#'
#' @importFrom magrittr %>%
#' @importFrom igraph graph_from_adjacency_matrix as.undirected
#' @importFrom parallel mclapply detectCores
#' @export
ensemble_cgbn <- function(dat,
                          discrete_variable = NULL,
                          num_iteration = 1,
                          boot = FALSE,
                          sub_ratio = 0.9,
                          sample_class = NULL,
                          hugin = FALSE,
                          constraints = NULL,
                          alpha = 0.01,
                          tol = 1e-04,
                          maxit = 0,
                          n_cores = NULL) {

  ## RHugin only when used
  if (!requireNamespace("RHugin", quietly = TRUE)) {
    stop(
      "This function requires the optional package 'RHugin', which is not installed.\n",
      "Install it and retry. macOS instructions: ",
      "https://rhugin.r-forge.r-project.org/InstallingRHuginMacOSX.html",
      call. = FALSE
    )
  }

  ## Input checks
  if (!(is.data.frame(dat) || is.matrix(dat))) stop("'dat' must either be a data.frame or matrix")
  if (is.null(colnames(dat))) stop("'dat' must have column names")
  if (is.null(rownames(dat))) stop("'dat' must have row names")
  if (!is.null(discrete_variable) && !is.character(discrete_variable)) stop("'discrete_variable' must be a character vector")
  if (!is.null(discrete_variable) && !all(discrete_variable %in% colnames(dat))) stop("'discrete_variable' must be in the columns of `dat`")
  if (!is.null(sample_class) && !is.atomic(sample_class)) stop("'sample_class' must be an atomic vector")
  if (!isTRUE(boot)) {
    if (is.null(sub_ratio) || !is.numeric(sub_ratio) || sub_ratio <= 0 || sub_ratio > 1)
      stop("'sub_ratio' must be in (0,1] when boot = FALSE")
  }
  if (!is.null(sample_class) && length(sample_class) != nrow(dat))
    stop("'sample_class' length must match nrow(dat)")

  ## Prepare input sample_data
  if (!is.null(sample_class)) {
    sample_data <- data.frame(ID = rownames(dat), class = as.factor(sample_class))
  } else {
    sample_data <- data.frame(ID = rownames(dat))
  }

  ## One iteration
  run_iteration <- function(i) {
    bagging <- resample(sample_data = sample_data,
                        sample_class = sample_class,
                        boot = boot,
                        sub_ratio = sub_ratio)

    resampled_data <- dat[bagging, , drop = FALSE]

    cgnetwork <- .fit_cgbn(dat = as.data.frame(resampled_data),
                           discrete_variable = discrete_variable,
                           constraints = constraints,
                           alpha = alpha,
                           tol = tol,
                           maxit = maxit)

    ## robust Hugin graph -> adjacency matrix -> igraph
    adj_mat <- .hugin_to_adjmat(cgnetwork$network)
    ig_net  <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed")
    ig_skel <- igraph::as.undirected(ig_net, mode = "collapse")

    res <- list(
      bagging     = bagging,
      ig_network  = ig_net,
      ig_skeleton = ig_skel
    )
    if (isTRUE(hugin)) res$hugin_network <- cgnetwork
    return(res)
  }

  ## Cores
  det <- parallel::detectCores()
  if (is.na(det)) det <- 2
  cores <- if (is.null(n_cores)) max(1, det - 1) else max(1, as.integer(n_cores))

  ## Parallel (Unix/macOS). On Windows, mclapply runs sequentially.
  results <- parallel::mclapply(seq_len(num_iteration), run_iteration, mc.cores = cores)

  ## Collect (named lists: iter_1, iter_2, ...)
  ig_networks  <- list()
  ig_skeletons <- list()
  baggings     <- list()
  if (isTRUE(hugin)) hugin_networks <- list()

  for (i in seq_along(results)) {
    nm <- paste0("iter_", i)
    baggings[[nm]]     <- results[[i]]$bagging
    ig_networks[[nm]]  <- results[[i]]$ig_network
    ig_skeletons[[nm]] <- results[[i]]$ig_skeleton
    if (isTRUE(hugin)) hugin_networks[[nm]] <- results[[i]]$hugin_network
  }

  if (isTRUE(hugin)) {
    return(list(
      ig_networks    = ig_networks,
      ig_skeletons   = ig_skeletons,
      baggings       = baggings,
      hugin_networks = hugin_networks
    ))
  } else {
    return(list(
      ig_networks   = ig_networks,
      ig_skeletons  = ig_skeletons,
      baggings      = baggings
    ))
  }
}

#' Learn a conditional Gaussian Bayesian network using PC-algorithm
#'
#' @keywords internal
#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate_if filter pull rename arrange
#' @importFrom tibble column_to_rownames
.fit_cgbn <- function(dat,
                      discrete_variable = NULL,
                      constraints = NULL,
                      alpha = 0.01,
                      tol = 1e-04,
                      maxit = 0) {

  if (!requireNamespace("RHugin", quietly = TRUE)) {
    stop("'.fit_cgbn()' requires 'RHugin'. See: https://rhugin.r-forge.r-project.org/InstallingRHuginMacOSX.html",
         call. = FALSE)
  }

  ## checks
  if (!is.data.frame(dat)) stop("'dat' must be a data.frame")
  if (is.null(colnames(dat))) stop("'dat' must have column names")
  if (!is.null(discrete_variable) && !is.character(discrete_variable)) stop("'discrete_variable' must be a character vector")
  if (!is.null(discrete_variable) && !all(discrete_variable %in% colnames(dat))) stop("'discrete_variable' must be in the columns of `dat`")
  if (is.null(rownames(dat))) stop("'dat' must have row names")

  if (!is.null(discrete_variable)) {
    discrete   <- dat %>% dplyr::select(discrete_variable)
    continuous <- dat %>% dplyr::select(-discrete_variable)

    ## ensure classes
    if (!all(sapply(discrete, is.factor))) {
      discrete <- discrete %>% dplyr::mutate_if(function(x) !is.factor(x), as.factor)
    }
    if (!all(unname(sapply(continuous, is.numeric))))
      stop("column vectors of 'continuous' should be numeric")

    num_dis <- ncol(discrete)
    num_con <- ncol(continuous)
    Data <- cbind(discrete, continuous)
  } else {
    num_dis <- 0
    num_con <- ncol(dat)
    if (!all(unname(sapply(dat, is.numeric))))
      stop("column vectors of 'continuous' should be numeric")
    Data <- dat
  }

  ## RHugin domain
  network <- RHugin::hugin.domain()
  num_nodes <- ncol(Data)

  node_df <- data.frame(
    node   = colnames(Data),
    class  = unname(sapply(Data, class)),
    levels = unname(sapply(sapply(Data, levels), length)),
    type   = c(rep("discrete", num_dis), rep("continuous", num_con)),
    stringsAsFactors = FALSE
  )

  ## Add nodes
  for (i in seq_len(num_nodes)) {
    if (node_df[i, "type"] == "discrete") {
      RHugin::add.node(domain = network,
                       name   = node_df[i, "node"],
                       states = levels(Data[[ node_df[i, "node"] ]]))
    } else {
      RHugin::add.node(domain = network,
                       name = node_df[i, "node"],
                       kind = "continuous")
    }
  }

  if (length(RHugin::get.nodes(network)) != num_nodes)
    stop("Number of nodes in Hugin domain does NOT match the input data!")

  ## Add cases
  RHugin::set.cases(network, Data)

  ## Structure learning (PC)
  if (!is.null(constraints)) {
    RHugin::learn.structure(network, alpha = alpha, constraints = constraints)
  } else {
    RHugin::learn.structure(network, alpha = alpha)
  }

  ## Touch experience tables (required before CPT learning for discrete nodes)
  for (i in seq_len(num_nodes)) {
    RHugin::get.table(network, node_df[i, "node"], type = "experience")
  }

  ## Compile & learn CPTs (EM)
  RHugin::compile(network)
  RHugin::learn.cpt(network, tol = tol, maxit = maxit)

  ## Marginals
  if (!is.null(discrete_variable)) {
    discrete_marginal   <- .get_discrete_marginal(network, node_df)
    continuous_marginal <- .get_continuous_marginal(network, node_df)
    output <- list(network = network,
                   node_df = node_df,
                   discrete_marginal = discrete_marginal,
                   continuous_marginal = continuous_marginal)
  } else {
    continuous_marginal <- .get_continuous_marginal(network, node_df)
    output <- list(network = network,
                   node_df = node_df,
                   continuous_marginal = continuous_marginal)
  }
  return(output)
}

# --- helpers ---------------------------------------------------------------

#' Coerce an RHugin graph to an adjacency matrix
#' @importFrom igraph graph_from_graphnel as_adjacency_matrix
#' @noRd
#' @keywords internal
.hugin_to_adjmat <- function(network) {

  gobj <- RHugin::as.graph.RHuginDomain(network)

  # case 1: RHugin already returns a matrix
  if (is.matrix(gobj)) {
    return(gobj)
  }

  # case 2: RHugin returns a 'graphNEL'
  if (inherits(gobj, "graphNEL")) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' is required for adjacency conversion.", call. = FALSE)
    }
    ig <- igraph::graph_from_graphnel(gobj)
    adj <- igraph::as_adjacency_matrix(ig, sparse = FALSE)
    return(adj)
  }

  stop("Unrecognized object returned by RHugin::as.graph.RHuginDomain().")
}

#' Obtain continuous marginals
#' @keywords internal
#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull select rename mutate arrange
#' @importFrom tibble column_to_rownames
.get_continuous_marginal <- function(network, node_df) {
  continuous_nodes <- node_df %>%
    dplyr::filter(type == "continuous") %>%
    dplyr::pull(var = "node")

  marg_list <- lapply(continuous_nodes, function(nm) {
    unlist(RHugin::get.marginal(network, nm))
  })
  marginal <- do.call(rbind, marg_list)

  marginal <- as.data.frame(marginal) %>%
    dplyr::select(mean, cov) %>%
    dplyr::rename(var = cov) %>%
    dplyr::mutate(node = continuous_nodes) %>%
    dplyr::arrange(node) %>%
    tibble::column_to_rownames(var = "node")

  return(marginal)
}

#' Obtain discrete marginals
#' @keywords internal
#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
.get_discrete_marginal <- function(network, node_df) {
  discrete_nodes <- node_df %>%
    dplyr::filter(type == "discrete") %>%
    dplyr::pull(var = "node")

  maxLevels <- as.numeric(max(node_df[node_df$node %in% discrete_nodes, ]$levels))
  marginal <- matrix(nrow = length(discrete_nodes),
                     ncol = maxLevels,
                     dimnames = list(discrete_nodes, sprintf("State %d", 1:maxLevels)))

  for (i in seq_along(discrete_nodes)) {
    tmp <- RHugin::get.marginal(network, discrete_nodes[i])
    marginal[i, 1:length(tmp$table$Freq)] <- tmp$table$Freq
  }
  return(as.data.frame(marginal))
}

## (optional) silence NSE NOTES on CRAN for dplyr column names
utils::globalVariables(c("mean", "cov", "node", "type"))
