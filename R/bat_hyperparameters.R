#' Run Phylogenetic Gaussian Process Hyperparameter Estimation
bat_hyperparameter_estimation <- function(df = smooth_and_regularise_call_densities(), phylogenetic_tree = phylogeny, basis = get_independent_components(), samples_per_mean = 4, n_mean_estimates = 4, n_samples = 10){

  traits <- sample_mean_estimates(df, samples_per_mean, n_mean_estimates, n_samples)

  basis_weights <- apply(traits, c(2,3), linear_coefficients, predictors = basis)
  species_names <- colnames(basis_weights)
  colnames(basis_weights) <- strtrim(species_names, 4)

  full_tree <- expand_tree(phylogenetic_tree, n_nodes = n_mean_estimates)

  hyperparameters <- apply(basis_weights[, full_tree$tip.label, ], c(1,3) , pou_type2mle, phylogenetic_tree = full_tree, optim_function = 'uobyqa', logl_function = pou_logl_fast, upper_initialisation = c(2, 50, 2), n_restarts = 5)

  return(hyperparameters)
}


#' Return coefficients for linear model
linear_coefficients <- function(response, predictors){
  mod <- lm(response ~ predictors)

  return(mod$coefficients)
}
