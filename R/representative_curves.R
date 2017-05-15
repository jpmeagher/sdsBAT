#' Representative spectral curves
#'
#' @export
species_mean <- function(df = smooth_and_regularise_call_densities(), ignored_columns = c(1,2), vector_length = length(df$chirp1[[1]])){
  individual_representation <- individual_means(df = df, ignored_columns = ignored_columns, vector_length = vector_length)
  species_representation <- tapply(individual_representation, unlist(df$species), list_vector_mean)

  return(species_representation)
}

individual_means <- function(df = smooth_and_regularise_call_densities(), ignored_columns = c(1,2), vector_length = length(df$chirp1[[1]])){
  considered_cells <- df[-ignored_columns]
  individual_means <- apply(considered_cells, 1, unlist)
  individual_means <- lapply(individual_means, na.omit)
  individual_means <- sapply(individual_means, call_mean, vector_length = vector_length)

  return(individual_means)
}


call_mean <- function(observations, vector_length = 208){
  observations <- matrix(observations, nrow = vector_length)
  call_mean <- apply(observations, 1, mean)

  return(list(call_mean))
}

list_vector_mean <- function(some_list, vector_length = 208){
  list_as_matrix <- matrix(unlist(some_list), nrow = vector_length)
  mean_vector <- apply(list_as_matrix, 1, mean)
  return(mean_vector)
}

