#' Calls per individual in Bat Echolocation call Dataset
#'
#' @param df Defaults to the Bat Echolocation call dataset
#'
#' @export
calls_per_individual <- function(df = calls){
  call_mask <- apply(df, c(1,2), function(x) !is.na(x) & length(x[[1]]) != 1)
  per_individual <- apply(call_mask, 1, sum)

  return(per_individual)
}

#' Sampling balanced by species
#'
#' @export
balanced_samples <- function(df = calls, categorical_variable = calls$species, individuals_per_category = 4, n_samples = 10, ignored_variables = c(1,2)){
  samples_per_individual <- calls_per_individual(df)
  categorical <- unlist(categorical_variable)
  categories <- unique(categorical)
  n_individuals <- nrow(df)

  samples <- matrix(rep(0), nrow = individuals_per_category*length(categories), ncol = n_samples)
  samples <- data.frame(samples)

  sampled_variables <- df[, -ignored_variables]

  for(j in 1:n_samples){
    individuals_index <- unlist(tapply(1:n_individuals, categorical, sample, individuals_per_category))
    variable_index <- sapply(samples_per_individual[individuals_index], sample.int, size = 1)
    index <- cbind(individuals_index, variable_index)

    balanced <- sampled_variables[index]
    samples[ , j] <- as.data.frame(I(balanced))
  }

  rownames(samples) <- names(individuals_index)

  return(samples)
}

#' Sample of Species Mean Estimates
#' @export
get_mean_estimates <- function(df = smooth_and_regularise_call_densities(), samples_per_mean = 4, n_mean_estimates = 4){
  temp_sample <- balanced_samples(df, individuals_per_category = n_mean_estimates, n_samples = samples_per_mean)
  full_sample <- array(unlist(temp_sample), dim = c(length(temp_sample[[1]][[1]]), dim(temp_sample)))
  colnames(full_sample) <- rownames(temp_sample)
  mean_estimate <- apply(full_sample, c(1,2), mean)

  return(mean_estimate)
}

#' Resampled estimates of Species mean
#' @export
sample_mean_estimates <- function(df = smooth_and_regularise_call_densities(), samples_per_mean = 4, n_mean_estimates = 4, n_samples = 10){
  temp_sample <- get_mean_estimates(df, samples_per_mean, n_mean_estimates)
  full_sample <- array(NA, dim = c(dim(temp_sample), n_samples))
  colnames(full_sample) <- colnames(temp_sample)
  for(i in 1:n_samples){
    full_sample[,,i] <- get_mean_estimates(df, samples_per_mean, n_mean_estimates)
  }

  return(full_sample)
}
