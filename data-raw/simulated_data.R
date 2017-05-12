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
  hyperparameter_estimates_simulations[, toString(i), ] <- suppressMessages(sdsBAT::pou_simulation_study(ln_hyperparameters = log(theta), phylogenetic_tree = tree, n_samples = N, logl_function = sdsBAT::pou_logl_fast, optim_function = 'uobyqa', lower_initialisation = c(0, 0, 0), upper_initialisation = c(2, 150, 2), n_restarts = 10)[1:3, ] )
  remove(tree)
}

#devtools::use_data(hyperparameter_estimates_simulations, overwrite = TRUE)

# Set Hyperparameter values
s_p <- 0.5
l <- 250
s_n <- 0.1

theta <- c(s_p, l, s_n)
# Set number of samples
N <- 1000

# Set number of observations at leaf nodes
n <- c(2^0, 2^1, 2^2, 2^3)

# Array of Hyperparameter estimates for simulated datasets
simulated_bat_estimation <- array(0, dim = c(length(theta), length(n), N))
rownames(simulated_bat_estimation) <- c('s_p', 'l', 's_n')
colnames(simulated_bat_estimation) <- n

# Set random seed
set.seed(3)

# Perform Simulation Study
for(i in n){
  tree <- sdsBAT::expand_tree(sdsBAT::phylogeny, n_nodes = i, evolutionary_time = 0)
  simulated_bat_estimation[, toString(i), ] <- suppressMessages(sdsBAT::pou_simulation_study(ln_hyperparameters = log(theta), phylogenetic_tree = tree, n_samples = N, logl_function = sdsBAT::pou_logl_fast, optim_function = 'uobyqa', lower_initialisation = c(0, 0, 0), upper_initialisation = c(2, 150, 2), n_restarts = 10)[1:3, ] )
  remove(tree)
}

#devtools::use_data(simulated_bat_estimation)

# Set Hyperparameter values
s_p <- 0.5
l <- 250
s_n <- 0.1

theta <- c(s_p, l, s_n)
# Set number of samples
N <- 1000

full_tree <- sdsBAT::expand_tree(sdsBAT::phylogeny, n_nodes = 50, evolutionary_time = 0)
resampling_tree <- sdsBAT::expand_tree(sdsBAT::phylogeny, n_nodes = 4, evolutionary_time = 0)

full_sample <- pou_simulated_data(log(theta), full_tree, 1)
plot(full_sample)
full_theta_est <- pou_type2mle(full_sample, full_tree, logl_function = pou_logl_fast, optim_function = 'uobyqa', upper_initialisation = c(2,300,2), n_restarts = 10 )
full_theta_est

resample_phylogeny <- function(full_sample, size = 4){
  resampled <- unlist(tapply(full_sample, rownames(full_sample), sample, size = size))
  return(resampled)
}

species_order <- cbind(resampling_tree$tip.label, rep(1:4, 22))
species_order <- apply(species_order, 1, paste, collapse = "")

resampled <- matrix(NA, nrow = 88, ncol = N)
for(i in 1:N){
  resampled[,i] <- resample_phylogeny(full_sample)[species_order]
}

resampled_hyperparameter_simulation <-suppressMessages(apply(resampled, 2, pou_type2mle, phylogenetic_tree = resampling_tree, logl_function = pou_logl_fast, optim_function = 'uobyqa', upper_initialisation = c(2,300,2), n_restarts = 10))

apply(resampled_hyperparameter_simulation, 1, mean)
apply(resampled_hyperparameter_simulation, 1, sd)
pairs(t(resampled_hyperparameter_simulation), diag.panel = panel.hist)

# devtools::use_data(full_theta_est)
# devtools::use_data(resampled_hyperparameter_simulation)
