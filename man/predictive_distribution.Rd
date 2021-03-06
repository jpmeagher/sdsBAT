% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bat_hyperparameters.R
\name{predictive_distribution}
\alias{predictive_distribution}
\title{Predictive Distribution for internal nodes}
\usage{
predictive_distribution(observed_weight = representative_basis_weights(),
  pou_hyperparameters = apply(bat_hyperparameter_estimates, c(1, 2), mean),
  phylogenetic_tree = phylogeny)
}
\description{
Predictive Distribution for internal nodes
}
