#' Perform community detection on a network
#'
#' @param ig an igraph object
#' @param mode the connection mode used to defined "isolated" nodes
#' @param vertex_symbol the node attribute that represents the symbol of the node
#' @param method the community detection method
#' @param small_community the threshold of small communities to be removed, is NULL, small communities will not be removed
#' @param steps the number of steps in walk-trap algorithm
#' @param resolution the resolution in Louvain and/or Leiden algorithm
#' @param n_iterations the number of iteration in Leiden algorithm
#'
#' @return a list object
#' @import methods utils

#' @export
community_detection <- function(ig,
                                mode=c('total','in','out'),
                                vertex_symbol='name',
                                method = c('walk_trap','louvain','leiden'),
                                small_community = NULL,
                                steps = 10,
                                resolution = 1,
                                n_iterations = 5){

  stopifnot(is(ig, "igraph"))
  if(is.null(vertex_symbol)){
    V(ig)$name <- as.character(seq_len(igraph::vcount(ig)))
    vertex_symbol <- "name"
  }

  if (!is.null(vertex_symbol) && !(vertex_symbol %in% igraph::vertex_attr_names(ig))) {
    stop("'vertex_symbol' NOT in `vertex_attr` of `ig`")
  }

  ## match arguments
  match.arg(mode)
  match.arg(method)

  ## remove isolated vertices
  igs <- remove_isolated(ig,
                         mode=mode,
                         vertex_symbol=vertex_symbol)

  ## community detection
  if(method == 'walk_trap'){igs <- .walktrap(igs,
                                             vertex_symbol=vertex_symbol,
                                             steps=steps)}

  if(method == 'louvain'){igs <- .louvain(igs,
                                          vertex_symbol=vertex_symbol,
                                          resolution=resolution)}

  if(method == 'leiden'){igs <- .leiden(igs,
                                        vertex_symbol=vertex_symbol,
                                        resolution=resolution,
                                        objective_function = "modularity",
                                        n_iterations = n_iterations)}

  if(!is.null(small_community)){
    igs <- remove_small_community(igs, threshold = small_community)
  }

  ## colorize communities
  igs$ig <- colorize_community(igs$ig, community_attr_name = "community")

  return(igs)
}


#' Perform walk-trap community detection algorithm on a network
#'
#' @param igs a list containing an igraph object named "ig"
#' @param vertex_symbol the node attribute that represents the symbol of the node
#' @param steps the number of steps in walk-trap algorithm
#'
#' @return a list object

#' @import methods utils
#' @importFrom igraph is_weighted cluster_walktrap vertex_attr modularity V E
#' @importFrom clustAnalytics conductance

#' @keywords internal
.walktrap <- function(igs,
                      vertex_symbol='name',
                      steps = 10){

  set.seed(1)
  if(igraph::is_weighted(igs$ig)){
    weights <- igraph::E(igs$ig)$weights
  }else{
    weights <- NULL
  }
  ig_c <- igraph::cluster_walktrap(graph = igs$ig,
                                   weights = weights,
                                   steps = steps)

  igraph::V(igs$ig)$community <- ig_c$membership

  community_partition <- list()
  for(i in 1:length(unique(igraph::V(igs$ig)$community))){
    community_partition[[paste0('C', i)]] <- igraph::vertex_attr(igs$ig, vertex_symbol, which(igraph::V(igs$ig)$community == i))
  }

  igs$community_partition <- community_partition

  igs$modularity <- round(igraph::modularity(ig_c, ig_c$membership), 2)

  igs$conductance <- round(clustAnalytics::conductance(igs$ig, ig_c$membership), 2)

  return(igs)
}

#' Perform Louvain community detection algorithm on a network
#'
#' @param igs a list containing an igraph object named "ig"
#' @param vertex_symbol the node attribute that represents the symbol of the node
#' @param resolution resolution in Louvain algorithm
#'
#' @return a list object

#' @import methods utils
#' @importFrom igraph is_weighted cluster_louvain vertex_attr modularity V E
#' @importFrom clustAnalytics conductance

#' @keywords internal
.louvain <- function(igs,
                     vertex_symbol='name',
                     resolution = 1){

  set.seed(1)
  if(igraph::is_weighted(igs$ig)){
    weights <- igraph::E(igs$ig)$weights
  }else{
    weights <- NULL
  }
  ig_c <- igraph::cluster_louvain(igs$ig,
                                  weights = weights,
                                  resolution=resolution)
  igraph::V(igs$ig)$community <- ig_c$membership

  community_partition <- list()
  for(i in 1:length(unique(igraph::V(igs$ig)$community))){
    community_partition[[paste0('C', i)]] <- igraph::vertex_attr(igs$ig, vertex_symbol, which(igraph::V(igs$ig)$community == i))
  }

  igs$community_partition <- community_partition

  igs$modularity <- round(igraph::modularity(ig_c, ig_c$membership), 2)

  igs$conductance <- round(clustAnalytics::conductance(igs$ig, ig_c$membership), 2)

  return(igs)
}

#' Perform Leiden community detection algorithm on a network
#'
#' @param igs a list containing an igraph object named "ig"
#' @param vertex_symbol the node attribute that represents the symbol of the node
#' @param resolution resolution in Leiden algorithm
#' @param objective_function the objective function in Leiden algorithm
#' @param n_iterations number of iteration in Leiden algorithm
#'
#' @return a list object

#' @import methods utils
#' @importFrom igraph is_weighted cluster_leiden vertex_attr modularity V E
#' @importFrom clustAnalytics conductance

#' @keywords internal
.leiden <- function(igs,
                    vertex_symbol='name',
                    resolution = 1,
                    objective_function = "modularity",
                    n_iterations = 5){

  set.seed(1)
  if(igraph::is_weighted(igs$ig)){
    weights <- igraph::E(igs$ig)$weights
  }else{
    weights <- NULL
  }
  ig_c <- igraph::cluster_leiden(igs$ig,
                                 weights = weights,
                                 objective_function = objective_function,
                                 resolution_parameter = resolution,
                                 n_iterations = n_iterations)

  igraph::V(igs$ig)$community <- ig_c$membership

  community_partition <- list()
  for(i in 1:length(unique(igraph::V(igs$ig)$community))){
    community_partition[[paste0('C', i)]] <- igraph::vertex_attr(igs$ig, vertex_symbol, which(igraph::V(igs$ig)$community == i))
  }

  igs$community_partition <- community_partition

  igs$modularity <- round(ig_c$quality,2)

  igs$conductance <- round(clustAnalytics::conductance(igs$ig, ig_c$membership), 2)


  return(igs)

}

