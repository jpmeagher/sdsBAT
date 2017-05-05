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
