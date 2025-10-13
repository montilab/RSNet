#' Simulate network via erdos-renyi model
#'
#' @param p The number of vertices in the graph
#' @param e_prob A probability for drawing an edge between two arbitrary vertices
#' @param DAG If TRUE, simulate an directed acyclic graph, or DAG

#'
#' @return An igraph object
#'
#' @importFrom igraph sample_gnp is_dag
#'
#' @export
model_er <- function(p=300, e_prob=0.01, DAG=FALSE) {

  if(!DAG){
    g <- igraph::sample_gnp(n=p,
                            p=e_prob,
                            directed=FALSE)
  }else{
    g <- igraph::sample_gnp(n=p,
                            p=e_prob,
                            directed=TRUE)
    while(!igraph::is_dag(g)){
      g <- igraph::sample_gnp(n=p,
                              p=e_prob,
                              directed=TRUE)
    }
  }
  return(g)
}

#' Simulate network via small world model
#'
#' @param p The number of vertices in the graph
#' @param dim The dimension of the starting lattice
#' @param nei The neighborhood within which the vertices of the lattice will be connected
#' @param rw_prob The probability for rewiring

#'
#' @return An igraph object
#'
#' @importFrom igraph sample_smallworld
#'
#' @export
model_sw <- function(p=300, dim=1, nei=3, rw_prob=0.05) {
  g <- igraph::sample_smallworld(size=p,
                                 dim=dim,
                                 nei=nei,
                                 p=rw_prob)
  return(g)
}


#' Simulate network via preferential attachment model
#'
#' @param p The number of vertices in the graph
#' @param power The power of the preferential attachment
#' @param z_appeal The attractiveness of the vertices with no adjacent edges
#' @param DAG If TRUE, simulate an directed acyclic graph, or DAG

#'
#' @return An igraph object
#'
#' @importFrom igraph sample_pa is_dag
#'
#' @export
model_pa <- function(p=300, power=1, z_appeal=1, DAG=FALSE) {
  if(!DAG){
    g <- igraph::sample_pa(n=p,
                           power=power,
                           zero.appeal=z_appeal,
                           directed=FALSE)
  }else{
    ## try "repeat" and "break" if "DAG"
    g <- igraph::sample_pa(n=p,
                           power=power,
                           zero.appeal=z_appeal,
                           directed=TRUE)
    while(!igraph::is_dag(g)){
      g <- igraph::sample_pa(n=p,
                             power=power,
                             zero.appeal=z_appeal,
                             directed=TRUE)
    }
  }
  return(g)
}

#' Simulate network via modular preferential attachment model
#'
#' @param p The number of vertices in the graph
#' @param m The number of modules
#' @param q_hub The quantile for degree cutoff for defining modular hubs
#' @param m_links The number of links between modules
#' @param power The power of the preferential attachment
#' @param z_appeal The attractiveness of the vertices with no adjacent edges

#'
#' @return An igraph object
#'
#' @importFrom igraph V E edges degree sample_pa
#'
#' @export
model_mpa <- function(p=300, m=6, q_hub=0.95, m_links=6, power=1.7, z_appeal=1) {

  # Calculate module sizes
  #m_size <- as.integer( (p - m_links) / m)
  m_size <- as.integer(p/m)

  # Create scale-free modules
  graphs <- lapply(seq_len(m), function(x) {
    g_x <- igraph::sample_pa(m_size,
                             power=power,
                             zero.appeal=z_appeal,
                             directed=FALSE)

    igraph::V(g_x)$module <- x
    igraph::V(g_x)$hub <- igraph::degree(g_x) >= quantile(igraph::degree(g_x), q_hub)
    igraph::V(g_x)$target <- !V(g_x)$hub

    return(g_x)
  })

  # Merge modules into a single graph
  g <- Reduce("+", graphs)

  pairs <- sapply(seq_len(m_links), function(x) {

    # Randomly choose two different modules to connect
    random.modules <- sample(unique(igraph::V(g)$module), size=2, replace=F)

    m1 <- as.numeric(random.modules[1])
    m2 <- as.numeric(random.modules[2])

    # Identify hub node of module 1
    m1_hubs <- as.character(igraph::V(g)[igraph::V(g)$hub & igraph::V(g)$module == m1])
    m2_targets <- as.character(igraph::V(g)[igraph::V(g)$target & igraph::V(g)$module == m2])

    c(sample(m1_hubs, size=1), sample(m2_targets, size=1))
  })

  # Add in new links colored red
  igraph::E(g)$color <- "black"
  e <- igraph::edges(pairs)
  e$color <- "red"

  g <- g + e
  return(g)
}

#' Simulate network via stochastic block model
#'
#' @param p The number of vertices in the graph
#' @param m The number of modules
#' @param sizes Size of each module, default = equal-sized modules
#' @param p_in within-module edge density
#' @param p_out between-module edges

#'
#' @return An igraph object
#'
#' @importFrom igraph sample_sbm
#'
#' @export
model_sbm <- function(p=300, m=6, sizes = NULL, p_in=0.8, p_out=0.05) {

  if(is.null(sizes)){
    sizes <- rep(p/m, m)
  }
  if(!is.null(sizes)){
    if( !(length(sizes) == m) ) stop("'sizes' must be a vector of length equals the number of modules, m")
    if( !(sum(sizes) == p ) ) stop(" 'sizes' must be a vector that sum up to the number of nodes, p")
  }
  print(sizes)
  # Create the Stochastic Block Model graph
  ## Initialize p_out for all
  block_matrix <- matrix(p_out, nrow = m, ncol = m)

  ## Set p_in for intra-module connections
  diag(block_matrix) <- p_in

  g <- igraph::sample_sbm(p, block_matrix, sizes)

  return(g)
}


#' Simulate multivariate Gaussian data for an undirected graph model
#'
#' @param n The number of samples to generate data for
#' @param ig An igraph object, must be an undirected graph
#' @param seed A number to seed the random data

#' @importFrom BDgraph bdgraph.sim
#' @importFrom igraph is_directed as_adjacency_matrix

#' @return A data frame of multivariate Gaussian data
#' @export
model_sim <- function(n, ig, seed=1){

  if(igraph::is_directed(ig)) stop("'ig' must be an undirected graph")

  # Extract adjacency matrix
  adj <- igraph::as_adjacency_matrix(ig, sparse = FALSE)

  # Simulate multivariate gaussian data
  set.seed(seed)

  bdg <- BDgraph::bdgraph.sim(n=n, graph=adj, type="Gaussian", vis=F)


  data <- as.data.frame(bdg$data)

  row.names(data) <- sprintf("S%d",1:nrow(data))

  colnames(data) <- sprintf("X%d",1:ncol(data))
  return(data)
}

#' Simulate family-based correlated data
#'
#' @param num_families number of family
#' @param members_per_family number of individuals in each family
#' @param h_sq inheritability
#' @param Sigma the covariance matrix of variables
#'
#' @return an ExpressionSet object

#' @import methods utils
#' @importFrom stats rnorm
#' @importFrom magrittr %>% set_rownames
#' @importFrom kinship2 makekinship
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame

#' @export
sim_fam_data <- function(num_families=100,
                         members_per_family = 5,
                         h_sq = 0.5,
                         Sigma){

  ## check the input variance-covariance matrix
  stopifnot(is.matrix(Sigma))
  stopifnot(isSymmetric(Sigma))

  ## number of variables
  p <- nrow(Sigma)

  # Total number of individuals
  N <- num_families * members_per_family

  # Family and individual IDs
  famid <- rep(1:num_families, each = members_per_family)
  id <- 1:N

  # Assign parents dynamically
  father <- rep(0, N)
  mother <- rep(0, N)

  if (members_per_family > 2) {
    for (f in 1:num_families) {
      father_idx <- (f - 1) * members_per_family + 1  # First person = father
      mother_idx <- (f - 1) * members_per_family + 2  # Second person = mother
      offspring_idx <- seq((f - 1) * members_per_family + 3, f * members_per_family)  # Remaining are offspring

      father[offspring_idx] <- father_idx
      mother[offspring_idx] <- mother_idx
    }
  }

  # Compute kinship matrix, and relatedness matrix
  kmat <- kinship2::makekinship(famid, id, father, mother)
  Phi <- 2*kmat

  # pedigree data
  data_ped <- data.frame(fam_id = famid,
                         subject_id = id,
                         father_id = father,
                         mother_id = mother) %>%
    magrittr::set_rownames(paste("S",sprintf("%d",1:N), sep = ""))


  ## data simulation
  h_mat <- h_sq * diag(p)
  i_mat <- diag(N)

  ## Np x Np
  X_vec_mat <- kronecker(i_mat, Sigma) + kronecker((Phi-i_mat), h_mat%*%Sigma)
  eigen_obj <- eigen(X_vec_mat)
  lambda <- eigen_obj$values
  eigenv <- eigen_obj$vectors
  temp_m <- eigenv %*% diag(lambda^(1/2))

  ## N x p
  X <- t(matrix(temp_m  %*% rnorm(nrow(X_vec_mat), 0,1), nrow = nrow(Sigma)))
  rm(X_vec_mat)
  rm(eigen_obj)
  rm(lambda)
  rm(eigenv)
  rm(temp_m)
  row.names(X) <- paste("S",sprintf("%d",1:N), sep = "")
  colnames(X) <- paste("X",sprintf("%d",1:p), sep = "")

  ## Use expression set
  return(Biobase::ExpressionSet(assayData = t(X),
                                phenoData = Biobase::AnnotatedDataFrame(data_ped)))
}



