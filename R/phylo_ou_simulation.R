#' Phylogenetic Ornstein-Uhlenbeck Process Covariance
#'
#' Derives the Covariance matrix for points in a univariate Phylogenetic
#' Gaussian Process given the Ornsetin-Uhelenbeck process kernel. \eqn{k(x, x')
#' = \sigma_p^2 \exp(\frac{|x - x'|}{\ell})\) + \sigma_n^2}
#'
#' @param hyperparameters Vector of hyperparameters for an Ornstein-Uhlenbeck
#'   kernel in form \eqn{(\sigma_p, \ell, \sigma_n)}
#' @param phylogenetic_tree Phylogenetic tree of relationships between
#'   observations
#'
#' @return Outputs the covariance between obsevations taken at the tips of the
#'   phylogenetic tree.
#'
#' @examples  K <- phylo_ou_covar(c(1,20,1), phylogeny)
phylo_ou_cov <- function(hyperparameters, phylogenetic_tree){
  if(length(hyperparameters) != 3){
    stop("Hyperparameters for phylogenetic noise, phylogenetic length-scale and non-phylogenetic noise must be included.")
  }
  if(any(theta < 0)){
    stop("Hyperparameters must be greater than 0")
  }
  if(attributes(phylogenetic_tree)$class != 'phylo'){
    stop("phylogenetic tree must be an object of class 'phylo'.")
  }

  distaces <- ape::cophenetic.phylo(phylogentetic_tree)
  ntips <- length(phylogenetic_tree$tip.label)

  sp <- hyperparameters[1]^2
  l <- hyperparameters[2]
  sn <- hyperparameters[3]^2

  K <- matrix( rep(0, ntips^2) , ncol= ntips)
  tempK <- K
  for ( i in 1:ntips){
    for ( j in i:ntips){
      rr <- abs(X[i,j])
      K[i,j] <- s_f2 * exp(-rr/l)
    }
  }
  tempK <- K + t(K)
  K <- tempK + diag(ntips)*(sn + .0000001 - sp)

  return(K)
}

#' Simulated Realisations of a Phylogenetic Ornstein-Uhlenbeck Process
#'
#' Title says it all really.
#'
#' @inheritParams phylo_ou_cov
#' @param nsamples The number of samples to be simulated.
#'
#' @return Outputs a matrix of size \eqn{ntips x nsamples} containing
#'   \eqn{nsamples} somulations from a given Phylogenetic Ornstein-Uhlenbeck
#'   Gaussian Process.
phylo_ou_simulated_obs <- function(hyperparameters, phylogenetic_tree, nsamples){
  if(length(hyperparameters) != 3){
    stop("Hyperparameters for phylogenetic noise, phylogenetic length-scale and non-phylogenetic noise must be included.")
  }
  if(any(theta < 0)){
    stop("Hyperparameters must be greater than 0")
  }
  if(attributes(phylogenetic_tree)$class != 'phylo'){
    stop("phylogenetic tree must be an object of class 'phylo'.")
  }
  K <- phylo_ou_cov(hyperparameters, phylogenetic_tree)
  ntips <- length(phylogenetic_tree$tip.label)
  C <- chol(K)

  simulated_data <- C%*%matrix(rnorm(ntips*nsamples, 0, 1), nrow = ntips, ncol = nsamples)

  return(simulated_data)
}
