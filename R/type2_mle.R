#' Hyperparameter Estimation
#'
#' Type II Maximum Likelihood Estimation for hyperparameters of a Phylogenetic
#' Ornstein-Uhlenbeck Process.
#'
#' @inheritParams pou_covariance
#' @inheritParams pou_logl_slow
#' @param logl_function function calculating the negative log likelihood of the
#'   observations to be minimised. Either 'pou_logl_slow' or 'pou_logl_fast'.
#' @param optim_function String indicating which optimising function to use,
#'   'optim' or 'uobyqa' from the minqa package.
#' @param optim_method Method to be used with the optim optimiser.
#' @param lower_initialisation Lower bound for the uniform random variable used
#'   to generate initial values.
#' @param upper_initialisation Upper bound for the uniform random variable used
#'   to generate initial values.
#' @param n_restarts Number of different initialisation values to test.
#'
#' @return Vector of optimised hyperparameters and the corresponding negative
#'   log likelihood, in order (phylogenetic noise, phylogenetic length-scale,
#'   non-phylogenetic noise, negative log likelihood)
pou_type2mle <- function(phylogenetic_tree, observations, logl_function = pou_logl_slow, optim_function = 'optim', optim_method = 'CG', lower_initialisation = c(0,0,0), upper_initialisation = c(1,1,1), n_restarts = 1){

  output <- rep(Inf, 4)
  names(output) <= c('sp', 'l', 'sn', 'ml')

  for(i in 1:n_restarts){
    initial_values <- c()
    for(j in 1:3){
      initial_values[j] <- runif(1, lower_initialisation[j], upper_initialisation[j])
    }
    initial_values <- log(initial_values)

    if(optim_function == 'optim'){
      sol <- optim(initial_values, logl_function, phylogenetic_tree = phylogenetic_tree, observations = observations, method = optim_method)

      hyperparameters <- exp(sol$par)
      ml <- sol$val

      if(ml < output[4]) {
        output[1:3] <- hyperparameters
        output[4] <- ml
      }

    }else if(optim_function == 'uobyqa'){
      sol <- minqa::uobyqa(initial_values, logl_function, phylogenetic_tree = phylogenetic_tree, observations = observations)

      hyperparameters <- exp(sol$par)
      ml <- sol$fval

      if(ml < output[4]) {
        output[1:3] <- hyperparameters
        output[4] <- ml
      }
    }else{stop("Please use either 'optim' or 'uobyqa'.")}

  }
  return(output)
}


