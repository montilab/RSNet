#' Find graphlet degree vector (GDV) of each vertex and graphlet correlation matrix of the given network
#'
#' @param ig an igraph object
#' @param level node-level of graphlet
#' @param redundant_orbit if TRUE, keep the redundant orbits, and remove otherwise
#' @param include_gcm if TRUE, compute the GCM
#' @param correlation the correlation measure between each orbit, "spearman" correlation is applied by default

#' @importFrom magrittr %>% set_rownames
#' @importFrom igraph as_edgelist vcount
#' @importFrom dplyr mutate_all mutate select
#' @importFrom orca count4 count5

#' @references Yaveroğlu, Ömer Nebil, et al. "Revealing the hidden language of complex networks."

#' @return a list of two matrics, gdvm = graphlet degree vector matrix, gcm = graphlet correlation matrix
#' @export
gdvm_gcm <- function(ig,
                     level = c("4","5"),
                     redundant_orbit=FALSE,
                     include_gcm = FALSE,
                     correlation = "spearman"){

  ## Default arguments
  level <- match.arg(level)

  ## number of vertices
  n <- igraph::vcount(ig)

  # Ensure "name" attribute exist
  if(!("name" %in% igraph::vertex_attr_names(ig))){igraph::V(ig)$name <- 1:n}

  edge_list <- igraph::as_edgelist(ig,names = FALSE) %>%
    as.data.frame() %>%
    dplyr::mutate_all(as.integer)

  ## count GDV of each vertex, including orbits
  if(level == "4"){ gdvm <- orca::count4(edge_list) }

  if(level == "5"){ gdvm <- orca::count5(edge_list) }

  ## fill missing rows with zeros to ensure the output has the correct number of rows
  if(nrow(gdvm) < n){
    fill_zero <- matrix(0, nrow = n-nrow(gdvm), ncol = ncol(gdvm))
    colnames(fill_zero) <- colnames(gdvm)

    gdvm <- rbind(gdvm, fill_zero)
  }

  stopifnot(nrow(gdvm) == n)

  ## remove redundant orbits
  if(!redundant_orbit){
    gdvm <- gdvm[, !colnames(gdvm) %in% paste0("O",c(3,13,12,14,16,21,20,23,38,28,26,44,47,69,17,72,71))]
  }

  ## Add row.names
  row.names(gdvm) <- igraph::V(ig)$name

  ## create the graphlet correlation matrix (GCM)
  ## Need to be aware of "NAs", since some orbits may have 0 counts across all vertices
  if(include_gcm){
    gcm <- cor(gdvm, method = correlation)

    return(list(gdvm = gdvm,
                gcm = gcm))
  }else{
    return(gdvm)
  }
}
