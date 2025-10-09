#' Function to count signed graphlets
#'
#' @param ig an igraph object
#' @param n_cores number of cores for parrallel computing
#' @param redundant if TRUE, include the redundant orbits
#' @param include_gcm if TRUE, compute GCM

#' @importFrom magrittr %>%
#' @importFrom igraph simplify edge_attr_names E vcount V delete_vertex_attr as_adj_list as_data_frame degree get.edge.ids graph are.connected
#' @importFrom dplyr select mutate arrange row_number
#' @importFrom parallel mclapply detectCores
#' @importFrom stats cor setNames

#' @references Das, Apratim. Efficient enumeration of small graphlets and orbits.

#' @return a list of two matrics, gdvm = graphlet degree vector matrix, gcm = graphlet correlation matrix
#' @export
signed_gdvm_gcm <- function(ig,
                            n_cores = NULL,
                            redundant=FALSE,
                            include_gcm=FALSE){


  ## ensure "sign" in edge attribute
  stopifnot("sign" %in% igraph::edge_attr_names(ig))

  ## ensure "sign" takes value of 1 and -1
  stopifnot(all(igraph::E(ig)$sign %in% c(1, -1)))

  ## number of vertices
  n <- igraph::vcount(ig)

  ## Ensure "name" attribute = 1,2,...,n
  if("name" %in% igraph::vertex_attr_names(ig)){
    igraph::V(ig)$old_name <- igraph::V(ig)$name
    ig <- igraph::delete_vertex_attr(ig, "name")
  }
  igraph::V(ig)$name <- 1:n


  ## pre-compute adjacency list
  adj_list <- igraph::as_adj_list(ig)

  ## pre-compute edge_list
  edge_list <- igraph::as_data_frame(ig, what="edges") %>%
    dplyr::select(1:2, sign) %>%
    dplyr::mutate(sign = as.numeric(sign)) %>%
    as.matrix()

  # Create a data frame to sort by degree and tie-break by node ID
  rank_df <- data.frame(
    node_id = seq(1, n),
    degree = igraph::degree(ig)
  ) %>%
    dplyr::arrange(degree, node_id) %>%
    dplyr::mutate(rank = dplyr::row_number())

  # Create named vector: names = node IDs, values = ranks
  node_rank <- stats::setNames(rank_df$rank, rank_df$node_id)

  ## target graphlet
  target <- igraph::graph(edges = c(1,2,2,3), directed = FALSE)

  ## Intialized "core"
  if(is.null(n_cores)){
    cores <- parallel::detectCores()
  }else{
    cores <- n_cores+1
  }

  #############################################################################
  ## Count triangle orbits
  GOM_tri_pl <- parallel::mclapply(1:length(adj_list), function(v, n) {
    ## Initialization
    local_GOM <- matrix(0, nrow = n, ncol = 6,
                        dimnames = list(seq(1, n), paste0("O", 9:14)))

    ## neighbors of v
    v_nei <- adj_list[[v]]
    if (length(v_nei) < 2) return(local_GOM)

    # permutated neighbors
    per_nei <- combn(v_nei, 2)
    connected_idx <- which(apply(per_nei, 2, function(x){igraph::get.edge.ids(ig, x)}) != 0)
    if (length(connected_idx) > 0) {
      connected_nei <- per_nei[, connected_idx, drop=FALSE]
      for (i in 1:ncol(connected_nei)) {
        u <- connected_nei[1, i]
        w <- connected_nei[2, i]
        if (node_rank[u] > node_rank[v] && node_rank[w] > node_rank[v]) {
          local_GOM <- .count_tri_orbits(v, u, w, edge_list, local_GOM)
        }
      }
    }

    return(local_GOM)

  }, n = n, mc.cores = cores - 1)

  GOM_tri <- Reduce("+", GOM_tri_pl)

  #############################################################################
  ## Count 2-star orbits
  GOM_2star_pl <- parallel::mclapply(1:length(adj_list), function(v, n, ig) {
    local_GOM <- matrix(0, nrow = n, ncol = 7,
                        dimnames = list(seq(1, n), paste0("O", 2:8)))

    v_nei <- adj_list[[v]]
    if (length(v_nei) < 2) return(local_GOM)

    per_nei <- combn(v_nei, 2, simplify = FALSE)
    for (pair in per_nei) {
      u <- pair[1]
      w <- pair[2]
      if (!igraph::are.connected(ig, u, w)) {
        e1 <- igraph::get.edge.ids(ig, c(v, u), error = FALSE)
        e2 <- igraph::get.edge.ids(ig, c(v, w), error = FALSE)
        if (e1 != 0 && e2 != 0) {
          edge_list_sub <- matrix(c(v, u, igraph::E(ig)$sign[e1],
                                    v, w, igraph::E(ig)$sign[e2]),
                                  ncol = 3, byrow = TRUE)
          colnames(edge_list_sub) <- c("from", "to", "sign")
          local_GOM <- .count_2star_orbits(edge_list_sub, local_GOM)
        }
      }
    }
    return(local_GOM)
  }, n = n, ig = ig, mc.cores = cores - 1)

  GOM_2star <- Reduce("+", GOM_2star_pl)

  ##############################################################################
  ## O0 and O1
  ig_pos <- igraph::subgraph.edges(ig, igraph::E(ig)[sign == 1], delete.vertices = FALSE)
  ig_neg <- igraph::subgraph.edges(ig, igraph::E(ig)[sign == -1], delete.vertices = FALSE)
  GOM_edge <- cbind(O0 = igraph::degree(ig_pos), O1 = igraph::degree(ig_neg))


  ## output
  gdvm <- cbind(GOM_edge, GOM_2star, GOM_tri)

  ## Add row.names
  row.names(gdvm) <- igraph::V(ig)$old_name

  ## redundant orbits
  if(!redundant){
    redundant_orbits <- c("O10", "O12", "O14")
    gdvm <- gdvm[, setdiff(colnames(gdvm), redundant_orbits)]
  }

  ## Need to be aware of "NAs", since some orbits may have 0 counts across all vertices
  if(include_gcm){

    gcm <- stats::cor(gdvm, method = "spearman")
    return(list(gdvm = gdvm,
                gcm = gcm))
  }else{
    return(gdvm)
  }

}


#' Function to count triangle orbits
#'
#' @param v a vertex id
#' @param u a vertex id
#' @param w a vertex id
#' @param edge_list a matrix with 3 columns, the third column = "sign"
#' @param GOM the graphlet orbit matrix
#'
#' @return the updated graphlet orbit matrix


#' @keywords internal
.count_tri_orbits <- function(v, u, w, edge_list, GOM){

  ## edge_list of the triangle graphlet
  edge_list_tri <- edge_list[edge_list[,1] %in% c(v,u,w) & edge_list[,2] %in% c(v,u,w), ]
  stopifnot(nrow(edge_list_tri) == 3)

  ## type of signed triangle
  type <- sum(edge_list_tri[,"sign"]) ## 3, 1, -1, -3

  ## vertices in the triangle
  tri_v <- unique(c(edge_list_tri[,1], edge_list_tri[,2]))
  stopifnot(length(tri_v) == 3)
  if(type == 3){
    ## '+++'
    GOM[tri_v, "O9"] <- GOM[tri_v, "O9"] + 1
  }else if(type == -3){
    ## '---'
    GOM[tri_v, "O14"] <- GOM[tri_v, "O14"] + 1
  }else if(type == 1){
    ## '++-'
    for(node in tri_v){
      s_sum <- sum(edge_list_tri[edge_list_tri[,1] %in% node | edge_list_tri[,2] %in% node, "sign"])
      if(s_sum == 2){
        GOM[node, "O10"] <- GOM[node, "O10"] + 1
      }else{
        GOM[node, "O11"] <- GOM[node, "O11"] + 1
      }
    }
  }else if(type == -1){
    ## '+--'
    for(node in tri_v){
      s_sum <- sum(edge_list_tri[edge_list_tri[,1] %in% node | edge_list_tri[,2] %in% node, "sign"])
      if(s_sum == 0){
        GOM[node, "O12"] <- GOM[node, "O12"] + 1
      }else{
        GOM[node, "O13"] <- GOM[node, "O13"] + 1
      }
    }
  }

  return(GOM)
}

#' Function to count 2-star orbits
#'
#' @param edge_list a matrix with 3 columns, the third column = "sign"
#' @param GOM the graphlet orbit matrix
#'
#' @return an updated GOM
#' @keywords internal
.count_2star_orbits <- function(edge_list, GOM){
  ## edge_list of the 2-star graphlet
  stopifnot(nrow(edge_list) == 2)

  # edge of each vertex in this graphlet
  edge_counts <- table(c(edge_list[,1],edge_list[,2]))
  ## center vertex
  v_d2 <- as.numeric(names(edge_counts[edge_counts==2]))
  stopifnot(length(v_d2) == 1)
  ## peripheral vertices
  v_d1 <- as.numeric(names(edge_counts[edge_counts==1]))
  stopifnot(length(v_d1) == 2)

  ## type of signed 2-star
  type <- sum(edge_list[,"sign"]) ## 2, 0, -2

  if(type == 2){
    ## '++'
    GOM[v_d2, "O3"] <- GOM[v_d2, "O3"] + 1
    GOM[v_d1, "O2"] <- GOM[v_d1, "O2"] + 1
  }else if(type == -2){
    ## '--'
    GOM[v_d2, "O8"] <- GOM[v_d2, "O8"] + 1
    GOM[v_d1, "O7"] <- GOM[v_d1, "O7"] + 1
  }else if(type == 0){
    ## "+-'
    GOM[v_d2, "O6"] <- GOM[v_d2, "O6"] + 1
    for(node in v_d1){
      s <- edge_list[edge_list[,1] %in% node | edge_list[,2] %in% node, "sign"]
      if(s > 0){
        GOM[node, "O4"] <- GOM[node, "O4"] + 1
      }else{
        GOM[node, "O5"] <- GOM[node, "O5"] + 1
      }
    }

  }

  return(GOM)
}

