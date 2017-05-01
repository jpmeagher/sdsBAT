#' Phylogenetic Covariance
#'
#' Derives the Covariance matrix for points in a univariate Phylogenetic
#' Gaussian Process given the Ornsetin-Uhelenbeck process kernel.
#'
#' @param ln_hyperparameters Vector of the natural logarithm of hyperparameters
#'   for an Ornstein-Uhlenbeck kernel in order \eqn{log(phylogenetic noise,
#'   phylogenetic length-scale, non-phylogenetic noise)}
#' @param phylogenetic_tree Phylogenetic tree of relationships between
#'   observations
#'
#' @return Outputs the covariance between obsevations at the tips of the
#'   phylogenetic tree.
#'
#' @examples  K <- pou_cov(log(c(1,20,1)), phylogeny)
pou_covariance <- function(ln_hyperparameters, phylogenetic_tree){
  if(length(ln_hyperparameters) != 3){
    stop("Hyperparameters for phylogenetic noise, phylogenetic length-scale and non-phylogenetic noise must be included.")
  }
  if(attributes(phylogenetic_tree)$class != 'phylo'){
    stop("phylogenetic tree must be an object of class 'phylo'.")
  }

  distances <- ape::cophenetic.phylo(phylogenetic_tree)
  ntips <- length(phylogenetic_tree$tip.label)

  sp <- exp(2*ln_hyperparameters[1])
  l <- exp(ln_hyperparameters[2])
  sn <- exp(2*ln_hyperparameters[3])

  K <- matrix( rep(0, ntips^2) , ncol= ntips)
  tempK <- K
  for ( i in 1:ntips){
    for ( j in i:ntips){
      rr <- abs(distances[i,j])
      K[i,j] <- sp * exp(-rr/l)
    }
  }
  tempK <- K + t(K)
  K <- tempK + diag(ntips)*(sn + .0000001 - sp)

  return(K)
}

#' Simulate Phylogenetic Data
#'
#' Simulated Realisations of a Phylogenetic Ornstein-Uhlenbeck Process
#'
#' @inheritParams phylo_ou_cov
#' @param nsamples The number of samples to be simulated.
#'
#' @return Outputs a matrix of size \eqn{ntips x nsamples} containing
#'   \eqn{nsamples} simulations from a given Phylogenetic Ornstein-Uhlenbeck
#'   Gaussian Process.
pou_simulated_data <- function(ln_hyperparameters, phylogenetic_tree, n_samples){

  K <- pou_covariance(ln_hyperparameters, phylogenetic_tree)
  ntips <- length(phylogenetic_tree$tip.label)
  C <- chol(K)

  simulated_data <- C%*%matrix(rnorm(ntips*n_samples, 0, 1), nrow = ntips, ncol = n_samples)

  return(simulated_data)
}
