#' Find the graphlet degree distance between two vertices
#'
#' @param v_gdv the graphlet degree vector of the first node, denoted as v
#' @param u_gdv the graphlet degree vector of the second node, denoted as u
#' @param w the weight vector assigned to each orbit, w = 1 for all orbits by default
#' @param similarity if True, return the similarity score, and distance otherwise

#' @return the distance or similarity score

#' @references Milenković, Tijana, and Nataša Pržulj. "Uncovering biological network function via graphlet degree signatures."

#' @export
gdv_distance <- function(v_gdv, u_gdv, w=1, similarity=FALSE){

  ## check the length
  stopifnot(length(v_gdv) == length(u_gdv))

  ## check the weight vector
  if(length(w) > 1 && (length(w) != length(v_gdv))){stop("Invalid weight vector!")}
  if(length(w) == 1){w <- rep(w, length(v_gdv))}

  ## distance of each orbit
  numerator <- abs(log10(v_gdv+1) - log10(u_gdv+1))
  denominator <- log10(pmax(v_gdv, u_gdv) + 2)
  distances <- numerator/denominator

  if(!similarity){
    return(sum(distances)/sum(w))
  }else{
    return(1 -sum(distances)/sum(w) )
  }

}
