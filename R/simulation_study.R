#' Phylogenetic Simulation Study
#'
#' Type II Maximum Likelihood Estimation of hyperparameters for n_samples
#' realisations of a Phylogenetic Ornstein-Uhlenbeck Process with given
#' hyperparameters
#'
#' @inheritParams pou_covariance
#' @inheritParams pou_simulatied_data
#' @inheritParams pou_type2mle
#'
#' @return A \eqn{4 x n_samples} matrix. Hyperparameters are ordered in the rows
#'   as (phylogenetic noise, phylogenetic length-scale, non-phylogenetic noise).
#'   Each column is a separate sample.
pou_simulation_study <- function(ln_hyperparameters, phylogenetic_tree, n_samples, logl_function = pou_logl_slow, optim_function = 'optim', optim_method = 'CG', lower_initialisation = c(0,0,0), upper_initialisation = c(1,1,1), n_restarts = 1){

  S <- pou_simulated_data(ln_hyperparameters, phylogenetic_tree, n_samples)

  MLE <- apply(S, 2, pou_type2mle, phylogenetic_tree = phylogenetic_tree, logl_function = logl_function, optim_function = optim_function, optim_method = optim_method, lower_initialisation = lower_initialisation, upper_initialisation = upper_initialisation, n_restarts = n_restarts)

  rownames(MLE) <- c('Phylo', 'L-S', 'Non-Phylo','-logL')
  return(MLE)
}
