#' Plot Consensus Network
#'
#' @description
#' Visualizes a consensus network using visNetwork, with optional query-based
#' subsetting and optional confidence-interval display for partial correlations.
#'
#' @param ig an igraph object
#' @param query character vector of node names/symbols to center the subgraph on
#' @param order neighborhood order around the query nodes
#' @param vertex_symbol vertex attribute to use as node label/symbol
#' @param edge_label edge attribute to show as label; default is "pcor"
#' @param edge_width optional edge attribute for edge thickness
#' @param CI_show logical; if TRUE, show CI **only** when edge_label = "pcor"
#'   and both "lower_quantile" and "upper_quantile" are available
#' @param main title
#' @param query_color color for queried nodes
#' @param neighbor_color color for neighbor nodes
#' @param all_color color for all nodes when no query is given
#' @param hover logical
#' @param stabilization logical
#' @param igraph_layout layout name to pass to visIgraphLayout
#' @param layout_seed random seed for layout
#'
#' @return list with p (visNetwork object) and data (nodes/edges)
#'
#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate case_when
#' @importFrom igraph vertex_attr_names edge_attr_names as_data_frame
#'   neighborhood induced_subgraph vcount
#' @importFrom visNetwork visNetwork toVisNetworkData visGroups visOptions
#'   visLegend visEdges visInteraction visPhysics visIgraphLayout visLayout
#'
#' @export
plot_cn <- function(ig,
                    query = NULL,
                    order = 1,
                    vertex_symbol = "name",
                    edge_label = "pcor",
                    edge_width = NULL,
                    CI_show = FALSE,
                    main = "Title",
                    query_color   = "#FF9999",
                    neighbor_color = "#56B4E9",
                    all_color      = "#56B4E9",
                    hover = TRUE,
                    stabilization = FALSE,
                    igraph_layout = "layout_with_fr",
                    layout_seed = 42) {

  stopifnot(methods::is(ig, "igraph"))

  # vertex attribute checks
  if (!(vertex_symbol %in% igraph::vertex_attr_names(ig))) {
    stop("`vertex_symbol` must be a vertex attribute of `ig`.")
  }
  if (!("name" %in% igraph::vertex_attr_names(ig))) {
    stop("Vertex attribute 'name' is required on `ig`.")
  }

  # edge attribute checks
  edge_attrs <- igraph::edge_attr_names(ig)
  if (!(edge_label %in% edge_attrs)) {
    stop(sprintf("edge_label = '%s' not found in edge attributes of ig.", edge_label))
  }
  if (!is.null(edge_width) && !(edge_width %in% edge_attrs)) {
    stop(sprintf("edge_width = '%s' not found in edge attributes of ig.", edge_width))
  }

  ## --- subset to query neighborhood if needed
  if (!is.null(query)) {
    query_ids <- igraph::as_data_frame(ig, what = "vertices") %>%
      dplyr::filter(name %in% query | !!as.name(vertex_symbol) %in% query) %>%
      row.names()

    all_vids <- igraph::neighborhood(ig, order = order, nodes = query_ids, mode = "all") %>%
      lapply(as.numeric) %>%
      unlist() %>%
      unique()

    ig <- igraph::induced_subgraph(ig, vids = all_vids)
  }

  if (igraph::vcount(ig) == 0) {
    message("The query graph is empty.")
    return(invisible(NULL))
  }

  vn_data <- visNetwork::toVisNetworkData(ig, idToLabel = FALSE)

  ## --- assign node groups
  if (!is.null(query)) {
    vn_data$nodes <- vn_data$nodes %>%
      dplyr::mutate(
        group = dplyr::case_when(
          (id %in% query) | (!!as.name(vertex_symbol) %in% query) ~ "query",
          TRUE ~ "neighbor"
        )
      )
  } else {
    vn_data$nodes <- vn_data$nodes %>%
      dplyr::mutate(group = "all")
  }

  ## --- edge labels: always build plain first
  vn_data$edges <- vn_data$edges %>%
    dplyr::mutate(
      label = as.character(!!as.name(edge_label)),
      dashes = dplyr::case_when(
        !!as.name(edge_label) < 0 ~ TRUE,
        TRUE ~ FALSE
      )
    )

  ## --- CI logic (strict)
  has_lower <- "lower_quantile" %in% colnames(vn_data$edges)
  has_upper <- "upper_quantile" %in% colnames(vn_data$edges)

  # CI is allowed ONLY if:
  # 1) CI_show = TRUE
  # 2) edge_label == "pcor"
  # 3) both CI columns exist
  can_show_ci <- isTRUE(CI_show) &&
    identical(edge_label, "pcor") &&
    has_lower && has_upper

  if (can_show_ci) {
    # overwrite the label we just made
    vn_data$edges <- vn_data$edges %>%
      dplyr::mutate(
        label = paste0(pcor, ",[", lower_quantile, ",", upper_quantile, "]")
      )
  } else {
    # CI was requested but not applicable → warn, keep plain label
    if (isTRUE(CI_show)) {
      if (!identical(edge_label, "pcor")) {
        warning("CI_show = TRUE but edge_label != 'pcor'; showing plain edge labels.")
      } else if (!(has_lower && has_upper)) {
        warning("CI_show = TRUE but 'lower_quantile' and/or 'upper_quantile' are missing; showing plain edge labels.")
      }
    }
  }

  ## --- edge width
  if (!is.null(edge_width)) {
    vn_data$edges <- vn_data$edges %>%
      dplyr::mutate(value = !!as.name(edge_width))
  }

  ## --- build visNetwork object
  p <- visNetwork::visNetwork(vn_data$nodes, vn_data$edges, main = main) %>%
    visNetwork::visGroups(groupname = "query",    color = query_color) %>%
    visNetwork::visGroups(groupname = "neighbor", color = neighbor_color) %>%
    visNetwork::visGroups(groupname = "all",      color = all_color) %>%
    visNetwork::visOptions(
      highlightNearest = TRUE,
      nodesIdSelection = TRUE,
      selectedBy = "group"
    ) %>%
    visNetwork::visLegend() %>%
    visNetwork::visEdges(smooth = FALSE) %>%
    visNetwork::visInteraction(hover = hover) %>%
    visNetwork::visPhysics(stabilization = stabilization) %>%
    visNetwork::visIgraphLayout(layout = igraph_layout) %>%
    visNetwork::visLayout(randomSeed = layout_seed)

  list(p = p, data = vn_data)
}
