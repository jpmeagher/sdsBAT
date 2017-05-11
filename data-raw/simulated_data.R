# Set Hyperparameter values
s_p <- 1
l <- 50
s_n <- 1

theta <- c(s_p, l, s_n)
# Set number of samples
N <- 1000

# Set number of observations at leaf nodes
n <- c(2^0, 2^1, 2^2, 2^3)

# Array of Hyperparameter estimates for simulated datasets
hyperparameter_estimates_simulations <- array(0, dim = c(length(theta), length(n), N))
rownames(hyperparameter_estimates_simulations) <- c('s_p', 'l', 's_n')
colnames(hyperparameter_estimates_simulations) <- n

# Set random seed
set.seed(3)

# Perform Simulation Study
for(i in n){
  tree <- sdsBAT::expand_tree(sdsBAT::phylogeny, n_nodes = i, evolutionary_time = 0)
  hyperparameter_estimates_simulations[, toString(i), ] <- suppressMessages(sdsBAT::pou_simulation_study(ln_hyperparameters = log(theta), phylogenetic_tree = tree, n_samples = N, logl_function = sdsBAT::pou_logl_fast, optim_function = 'uobyqa', lower_initialisation = c(0, 0, 0), upper_initialisation = c(2, 150, 2), n_restarts = 5)[1:3, ] )
  remove(tree)
}

# devtools::use_data(hyperparameter_estimates_simulations, overwrite = TRUE)


