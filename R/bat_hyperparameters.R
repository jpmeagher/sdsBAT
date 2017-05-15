#' Return coefficients for linear model
linear_coefficients <- function(response, predictors){
  mod <- lm(response ~ predictors)

  return(mod$coefficients)
}

#' Run Phylogenetic Gaussian Process Hyperparameter Estimation
#' @export
bat_hyperparameter_estimation <- function(df = smooth_and_regularise_call_densities(), phylogenetic_tree = phylogeny, basis = get_independent_components(), samples_per_mean = 4, n_mean_estimates = 4, n_samples = 10){

  traits <- sample_mean_estimates(df, samples_per_mean, n_mean_estimates, n_samples)
  traits <- apply(traits, c(2,3), '-', spectral_mean)

  basis_weights <- apply(traits, c(2,3), linear_coefficients, predictors = basis)

  full_tree <- expand_tree(phylogenetic_tree, n_nodes = n_mean_estimates)

  species_order <- cbind(full_tree$tip.label, rep(1:4, 22))
  species_order <- apply(species_order, 1, paste, collapse = "")

  hyperparameters <- apply(basis_weights[, species_order, ], c(1,3) , pou_type2mle, phylogenetic_tree = full_tree, optim_function = 'uobyqa', logl_function = pou_logl_fast, upper_initialisation = c(2, 50, 2), n_restarts = 5)

  return(hyperparameters)
}

#' Representative weights for each species
#' @export
representative_basis_weights <- function(representative_curves = species_mean(), basis = get_independent_components()){
  representative_curves <- lapply(representative_curves, '-', spectral_mean)
  basis_weights <- sapply(representative_curves, linear_coefficients, predictors = basis)
  return(basis_weights)
}

full_tree_covariance <- function(ln_hyperparameters, phylogenetic_tree){
  if(length(ln_hyperparameters) != 3){
    stop("Hyperparameters for phylogenetic noise, phylogenetic length-scale and non-phylogenetic noise must be included.")
  }
  if(attributes(phylogenetic_tree)$class != 'phylo'){
    stop("phylogenetic tree must be an object of class 'phylo'.")
  }

  distances <- ape::dist.nodes(phylogenetic_tree)
  n_nodes <- dim(distances)[1]

  sp <- exp(2*ln_hyperparameters[1])
  l <- exp(ln_hyperparameters[2])
  sn <- exp(2*ln_hyperparameters[3])

  K <- matrix( rep(0, n_nodes^2) , ncol= n_nodes)
  tempK <- K
  for ( i in 1:n_nodes){
    for ( j in i:n_nodes){
      rr <- abs(distances[i,j])
      K[i,j] <- sp * exp(-rr/l)
    }
  }
  tempK <- K + t(K)
  K <- tempK + diag(n_nodes)*(sn + .0000001 - sp)

  return(K)
}

#' Predictive Distribution for internal nodes
#' @export
predictive_distribution <- function(observed_weight = representative_basis_weights(), pou_hyperparameters = apply(bat_hyperparameter_estimates, c(1,2), mean), phylogenetic_tree = phylogeny){

  K <- array(dim = c(max(phylogenetic_tree$edge), max(phylogenetic_tree$edge), nrow(observed_weight)))

  predictive_mean <- matrix(nrow = phylogenetic_tree$Nnode, ncol = nrow(observed_weight))
  predictive_covariance <- array(dim = c(phylogenetic_tree$Nnode, phylogenetic_tree$Nnode, nrow(observed_weight)))

  tip_index <- seq_along(phylogenetic_tree$tip.label)

  for(i in seq_along(observed_weight[,1])){

    K[,,i] <- full_tree_covariance(log(pou_hyperparameters[1:3,i]), phylogenetic_tree = phylogenetic_tree)

    observed_K_inv <- solve(K[tip_index, tip_index, i])

    predictive_mean[,i] <- (K[-tip_index, tip_index, i] %*% observed_K_inv) %*% observed_weight[i,]

    predictive_covariance[,,i] <- K[-tip_index, -tip_index, i] - (K[-tip_index, tip_index, i] %*% observed_K_inv) %*% K[tip_index, -tip_index, i]
  }

  pred_distribution <- list(predictive_mean = predictive_mean, predictive_covariance = predictive_covariance)

  return(pred_distribution)

}

#' Ancestral Reconstruction of Spectral Curves
#' @export
ancestral_reconstruction <- function(predicted_weights = predictive_distribution()$predictive_mean, basis = get_independent_components()){
  basis <- cbind(rep(1, nrow(basis)), basis)
  ancestral_curves <- predicted_weights %*% t(basis)
  ancestral_curves <- apply(ancestral_curves, 1, '+', spectral_mean)

  return(ancestral_curves)
}

