#' Plot Consensus Network
#'
#' @param ig an igraph object
#' @param query query nodes for subgraph
#' @param order the order of neighbors of the queried nodes
#' @param vertex_symbol the node attribute that represents the symbol
#' @param edge_label the edge attributes that represents the edge label
#' @param edge_width the edge attributes that represents the edge width
#' @param CI_show if TRUE, plot the confidence interval
#' @param main title
#' @param query_color the color of queried nodes
#' @param neighbor_color the color of the neighbors of queries nodes
#' @param all_color the color for all nodes in the network
#' @param hover a logical argument
#' @param stabilization a logical argument
#' @param igraph_layout the layout being used
#' @param layout_seed the random seed of the layout

#' @import methods utils
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate case_when
#' @importFrom igraph vertex_attr_names as_data_frame neighborhood induced_subgraph vcount
#' @importFrom visNetwork visNetwork toVisNetworkData visGroups visOptions visLegend visEdges visInteraction visPhysics visIgraphLayout visLayout

#' @return a list containing a visNetwork object

#' @export
plot_cn <- function(ig,
                    query=NULL,
                    order = 1,
                    vertex_symbol = "symbol",
                    edge_label = NULL,
                    edge_width = NULL,
                    CI_show = FALSE,
                    main = "Title",
                    query_color="#FF9999",
                    neighbor_color= "#56B4E9",
                    all_color = "#56B4E9",
                    hover = TRUE,
                    stabilization = FALSE,
                    igraph_layout = "layout_with_fr",
                    layout_seed = 42){

  stopifnot(is(ig,"igraph"))
  if(!(vertex_symbol %in% igraph::vertex_attr_names(ig))) stop("`vertex_symbol` should be in `vertex_attr` of `ig`")

  if(!is.null(edge_label)){
    if(!(edge_label %in% igraph::edge_attr_names(ig))) stop(" invalid 'edge_label'")
  }

  if(!is.null(edge_width)){
    if(!(edge_width %in% igraph::edge_attr_names(ig))) stop(" invalid 'edge_width'")
  }

  if(!("name" %in% igraph::vertex_attr_names(ig))) stop(" 'name' should be in `vertex_attr` of `ig`")

  ## query if sub-graph
  if(!is.null(query)){

    ## obtain node ids of query nodes
    query_ids <- igraph::as_data_frame(ig, what='vertices') %>%
      dplyr::filter(name %in% query | !!as.name(vertex_symbol) %in% query) %>%
      row.names(.)

    ## obtain neighbors of query nodes
    all_vids <- lapply(igraph::neighborhood(ig, order=order, nodes= query_ids, mode='all'), as.numeric) %>%
      unlist(.) %>%
      unique(.)

    ## Obtain the node-induced sub-graph
    ig <- igraph::induced_subgraph(ig, vids= all_vids)
  }

  if(igraph::vcount(ig) == 0){
    print("The query graph is empty")
    return()
  }

  vn_data <- visNetwork::toVisNetworkData(ig, idToLabel=FALSE)


  # vn_data$nodes <- vn_data$nodes %>% dplyr::mutate(label = case_when(
  #   is.na(!!as.name(vertex_symbol)) ~ id,
  #   TRUE ~ !!as.name(vertex_symbol)))

  ## if query, group = "query" and "neighbor"
  if(!is.null(query)){
    vn_data$nodes <- vn_data$nodes %>%
      dplyr::mutate(group = case_when(
        (id %in% query)|(!!as.name(vertex_symbol) %in% query) ~ "query",
        TRUE ~ "neighbor"
      ))

  }

  ## plot the entire network, group = "all"
  if(is.null(query)){
    vn_data$nodes <- vn_data$nodes %>%
      dplyr::mutate(group = "all")
  }

  ## labeling edges with confidence intervals
  if(!is.null(edge_label) & CI_show ){
    if(!all(c("lower_quantile","upper_quantile") %in% colnames(vn_data$edges))) stop("Valid confidence interval NOT available, please label them as 'lower_quantile' and 'upper_qunatile' or set 'CI_show' as FALSE ")
    vn_data$edges <- vn_data$edges %>%
      dplyr::mutate(label = paste(!!as.name(edge_label), paste0("[",paste(lower_quantile,upper_quantile, sep=","),"]"),sep=",")) %>%
      dplyr::mutate(dashes = case_when(
        !!as.name(edge_label) < 0 ~ TRUE,
        TRUE ~ FALSE
      ))
  }
  if(!is.null(edge_label) & !CI_show ){
    vn_data$edges <- vn_data$edges %>%
      dplyr::mutate(label = as.character(!!as.name(edge_label))) %>%
      dplyr::mutate(dashes = case_when(
        !!as.name(edge_label) < 0 ~ TRUE,
        TRUE ~ FALSE
      ))
  }

  ## Assign edge thickness
  if(!is.null(edge_width)){
    vn_data$edges <- vn_data$edges %>%
      dplyr::mutate(value = !!as.name(edge_width))
  }

  p <- visNetwork::visNetwork(vn_data$nodes, vn_data$edges, main=main) %>%
    visNetwork::visGroups(groupname = "query", color = query_color) %>%
    visNetwork::visGroups(groupname = "neighbor", color = neighbor_color) %>%
    visNetwork::visGroups(groupname = "all", color = all_color) %>%
    visNetwork::visOptions(highlightNearest = TRUE,
                           nodesIdSelection = TRUE,
                           selectedBy = "group") %>%
    visNetwork::visLegend() %>%
    visEdges(smooth=FALSE) %>%
    visNetwork::visInteraction(hover = hover) %>%
    visNetwork::visPhysics(stabilization = stabilization) %>%
    visNetwork::visIgraphLayout(layout = igraph_layout) %>%
    visNetwork::visLayout(randomSeed = layout_seed)

  return(list(p = p,
              data = vn_data))
}
