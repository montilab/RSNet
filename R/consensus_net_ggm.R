#' Create consensus network of resampling-based Markov network
#'
#' @param ggm_networks the output from "ensemble_sggm.R"
#' @param CI confidence interval of interests
#' @param pos_cut a threshold for positive partial correlation
#' @param neg_cut a threshold for negative partial correlation
#' @param node_annot a data frame have at least two columns: 'id' and 'symbol', and 'id' must match the colnames in "mat"
#' @param filter filter method
#' @param threshold threshold of the selected filter


#' @import methods utils
#' @importFrom stats pnorm
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join pull
#' @importFrom tibble rownames_to_column
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist as_data_frame edge_attr E V

#' @return a list object including an igraph object as the consensus network
#' @export
consensus_net_ggm <- function(ggm_networks,
                              CI = 0.95,
                              pos_cut = 0,
                              neg_cut = 0,
                              node_annot = NULL,
                              filter = c("pval","fdr","none"),
                              threshold = 0.05){

  ## obtain method for statistical inference
  method <- ggm_networks$method
  p <- ggm_networks$p
  estimate_CI <- ggm_networks$estimate_CI

  # inference method
  match.arg(filter)
  if(!(filter %in% c("pval","fdr","none"))) stop("invalid filter!")

  if(!is.null(node_annot)){
    if(!is.data.frame(node_annot)) stop(" 'node_annot ' must be a data frame")
    if(!(c("id","symbol") %in% colnames(node_annot))) stop(" 'node_annot ' must contains 'id' and 'symbol' column ")
  }

  ## dealing with unexpected missing values
  ggm_networks$avg_partialCor[is.na(ggm_networks$avg_partialCor)] <- 0
  ggm_networks$avg_z_score_partialCor[is.na(ggm_networks$avg_z_score_partialCor)] <- 0

  ## estimate the p-value of each edge
  pval_pcor <- sapply(ggm_networks$avg_z_score_partialCor, .pvalue) %>%
    upper_tri_to_matrix(., variable_names = ggm_networks$vids, diagl = 1)

  ## estimate the adjusted p-value of each edge
  qval_pcor <- pval_pcor %>%
    matrix_p_adjust(.)

  ## Initialize consensus network
  if(!estimate_CI){
    pcor_avg <- ggm_networks$avg_partialCor %>%
      round(.,2) %>%
      upper_tri_to_matrix(., variable_names = ggm_networks$vids, diagl = 1)

    ## create consensus network, weight = mean partial correlation
    consensus_network <- apply(pcor_avg, c(1,2), .abs_pcor_filter, pos_cut, neg_cut)

    rm(ggm_networks)
  }else{
    pcor_avg <- apply(ggm_networks$partialCor_mat, 2, mean) %>%
      round(.,2) %>%
      upper_tri_to_matrix(., variable_names = ggm_networks$vids, diagl = 1)
    ## lower quantile
    CI_lower <- apply(ggm_networks$partialCor_mat, 2, quantile, (1-CI)/2) %>%
      round(.,2) %>%
      upper_tri_to_matrix(., variable_names = ggm_networks$vids, diagl = 1)
    ## upper quantile
    CI_upper <- apply(ggm_networks$partialCor_mat, 2, quantile, (1+CI)/2) %>%
      round(.,2) %>%
      upper_tri_to_matrix(., variable_names = ggm_networks$vids, diagl = 1)

    rm(ggm_networks)

    ## for export to excel format purpose
    pcor_CI <- paste(CI_lower, pcor_avg, CI_upper, sep = ",") %>%
      matrix(nrow = p, ncol = p)
    row.names(pcor_CI) <- row.names(pcor_avg)
    colnames(pcor_CI) <- colnames(pcor_avg)

    ## create consensus network, weight = mean partial correlation
    consensus_network <- apply(pcor_CI, c(1,2), .CI_filter)
    consensus_network <- apply(consensus_network, c(1,2), .abs_pcor_filter, pos_cut, neg_cut)
  }

  ## further filtering by significance
  if(filter == "fdr"){
    consensus_network[qval_pcor >= threshold] <- 0
  }
  if(filter == "pval"){
    consensus_network[pval_pcor >= threshold] <- 0
  }

  ## diagonal to 0 before creating igraph object
  diag(consensus_network) <- 0

  # to igraph object
  consensus_network  <- igraph::graph_from_adjacency_matrix(consensus_network,
                                                            mode = "undirected",
                                                            weighted = TRUE)

  names(igraph::edge_attr(consensus_network ))[which(names(igraph::edge_attr(consensus_network )) == "weight")] <- "pcor"
  igraph::E(consensus_network)$abs_pcor <- abs( igraph::E(consensus_network)$pcor)

  ## add extra node attributes
  if(!is.null(node_annot)){
    igraph::V(consensus_network)$symbol <- igraph::as_data_frame(consensus_network, what='vertices') %>%
      tibble::rownames_to_column(var='id') %>%
      dplyr::left_join(node_annot, by='id') %>%
      dplyr::pull(symbol)

    row.names(pcor_CI) <- igraph::V(consensus_network)$symbol
    colnames(pcor_CI) <- igraph::V(consensus_network)$symbol
  }

  if(estimate_CI){
    ## add extra edge attributes
    edgelist <- igraph::as_edgelist(consensus_network)

    igraph::E(consensus_network)$lower_quantile <- round(apply(edgelist, 1, .get_info, CI_lower),2)

    igraph::E(consensus_network)$upper_quantile <- round(apply(edgelist, 1, .get_info, CI_upper),2)

    igraph::E(consensus_network)$pval <- apply(edgelist, 1, .get_info, pval_pcor)

    igraph::E(consensus_network)$qval <- apply(edgelist, 1, .get_info, qval_pcor)

    igraph::E(consensus_network)$sign <- ifelse(E(consensus_network)$pcor > 0, 1, -1)

    return(list(consensus_network = consensus_network,
                pcor_CI = pcor_CI,
                CI = CI,
                method = method))
  }else{
    ## add extra edge attributes
    edgelist <- igraph::as_edgelist(consensus_network)

    igraph::E(consensus_network)$pval <- apply(edgelist, 1, .get_info, pval_pcor)

    igraph::E(consensus_network)$qval <- apply(edgelist, 1, .get_info, qval_pcor)

    igraph::E(consensus_network)$sign <- ifelse(E(consensus_network)$pcor > 0, 1, -1)

    return(list(consensus_network = consensus_network,
                method = method))
  }

}


#' Compute the two-sided p-value of a z-score
#'
#' @param z_score a z-score
#'
#' @return the p-value
#' @keywords internal
.pvalue <- function(z_score){

  return(2*stats::pnorm(q=abs(z_score), lower.tail=FALSE))
}


#' Filter edges whose effect size is small
#'
#' @param x the average partial correlation
#' @param pos_cut a threshold for positive partial correlation
#' @param neg_cut a threshold for negative partial correlation
#' @return the average value or 0
#' @keywords internal
.abs_pcor_filter <- function(x,
                             pos_cut,
                             neg_cut){

  if( (x) > 0 & (abs(x) <= abs(pos_cut)) ){
    return(0)
  }

  if( (x) < 0 & (abs(x) <= abs(neg_cut)) ){
    return(0)
  }

  return(x)
}

#' Convert the string of "low,avg,high" to numeric value
#'
#' @param x a string in the format of "low,avg,high"
#'
#' @return a numeric vector
#' @keywords internal

#' @importFrom magrittr %>%
#' @importFrom stringr str_split
.to_numeric=function(x){
  ## x = individual entry of the pcor_CI matrix
  x_to_numeric <- stringr::str_split(x,",") %>%
    unlist(.) %>%
    as.numeric(.)

  return(c(x_to_numeric[1], x_to_numeric[2], x_to_numeric[3] ))
}

#' Filter edges whose confidence interval contains 0
#'
#' @param x a string in the format of "low,avg,high"
#' @return the average value or 0
#' @keywords internal
.CI_filter=function(x){

  ## x = individual entry of the pcor_CI matrix
  x_to_numeric <- .to_numeric(x)

  ## lower quantile, mean, upper quantile
  low <- x_to_numeric[1]
  avg <- x_to_numeric[2]
  high <- x_to_numeric[3]


  ## check if confidence interval includes 0
  if(low > 0 & high > 0){
    avg <- avg
  }
  else if(low < 0 & high < 0){
    avg <- avg
  }
  else{return(0)}
}

#' Get the lower or upper quantile of an edge
#'
#' @param edge an edge of format (source, target)
#' @param info the information of interests
#'
#' @return the quantile information of the edge

#' @keywords internal
.get_info=function(edge, info){

  return(info[edge[1], edge[2]])
}
