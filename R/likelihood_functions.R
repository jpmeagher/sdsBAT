#' Phylogenetic Likelihood
#'
#' Calculates the Phylogenetic Ornstein-Uhelenbeck Process negative log
#' likelihood, given the process hyperparameters, a phylogenetic tree, and some
#' data.
#'
#' @param observations Observations from the tips of the phylogenetic tree
#' @inheritParams pou_covariance
#'
#' @return A positive value which is the negative log likelihood for the data.
#'
#' @export
pou_logl_slow <- function(ln_hyperparameters, phylogenetic_tree, observations){
  if (length(observations) != length(phylogenetic_tree$tip.label)){
    stop("Each observation must correspond to a single tip of the phylogenetic tree")
  }
  K <- pou_covariance(ln_hyperparameters, phylogenetic_tree)

  ntips <- length(phylogenetic_tree$tip.label)

  C <- chol(K)
  K_inv <- chol2inv(C)
  y <- matrix(observations, nrow = ntips, ncol = 1)

  A.1 <- -.5*((t(y) %*% K_inv) %*% y)
  A.2 <- -.5*log(prod(diag(C))^2)
  A.3 <- -(ntips/2) * log(2*pi)
  lik <- c(-(A.1 + A.2  + A.3))

  return(lik)
}

#' Fast Phylogenetic Likelihood
#'
#' Using the inline and RcppEigen packages performs a fast calculation of the
#' Phylogenetic Ornsetin-Uhelenbeck Process negative log likelihood, given the
#' process hyperparameters, a phylogenetic tree, and some data. The C++ code was
#' taken from code accompanying Hajipantelis et al.
#'
#' @inheritParams pou_logl_slow
#' @inheritParams pou_covariance
#'
#' @return A positive value which is the negative log likelihood for the data.
#'
#' @source P. Z. Hadjipantelis, N. S. Jones, J. Moriarty, D. A. Springate, and
#'   C. G. Knight, Function-valued traits in evolution, Journal of The Royal
#'   Society Interface. 10(82), 20121032 (2013). \url{https://github.com/fpgpr}
#'
#' @useDynLib sdsBAT
#' @importFrom Rcpp sourceCpp
#'
#' @export
pou_logl_fast <- function(ln_hyperparameters, phylogenetic_tree, observations){
  if(!requireNamespace("RcppEigen", quietly = TRUE)){
    stop("RcppEigen needed for this function to work. Please install it, or use pou_logl_slow()",
      call. = FALSE)
  }

  distances <- ape::cophenetic.phylo(phylogenetic_tree)

  storage.mode(distances) <- "double";        #Save stuff as doubles
  storage.mode(observations) <- "double";
  storage.mode(ln_hyperparameters) <- "double";

  return ( logl_cxx(ln_hyperparameters, observations, distances)$F);   #LogLikelihood
}


