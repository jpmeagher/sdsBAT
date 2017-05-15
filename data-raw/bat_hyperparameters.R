set.seed(3)

bat_hyperparameter_estimates <- bat_hyperparameter_estimation(n_samples = 1000)

devtools::use_data(bat_hyperparameter_estimates, overwrite = T)
