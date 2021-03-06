% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_study.R
\name{pou_simulation_study}
\alias{pou_simulation_study}
\title{Phylogenetic Simulation Study}
\usage{
pou_simulation_study(ln_hyperparameters, phylogenetic_tree, n_samples,
  logl_function = pou_logl_slow, optim_function = "optim",
  optim_method = "Nelder-Mead", lower_initialisation = c(0, 0, 0),
  upper_initialisation = c(1, 1, 1), n_restarts = 1)
}
\arguments{
\item{ln_hyperparameters}{Vector of the natural logarithm of hyperparameters
for an Ornstein-Uhlenbeck kernel in order \eqn{log(phylogenetic noise,
phylogenetic length-scale, non-phylogenetic noise)}}

\item{phylogenetic_tree}{Phylogenetic tree of relationships between
observations}

\item{logl_function}{function calculating the negative log likelihood of the
observations to be minimised. Either 'pou_logl_slow' or 'pou_logl_fast'.}

\item{optim_function}{String indicating which optimising function to use,
'optim' or 'uobyqa' from the minqa package.}

\item{optim_method}{Method to be used with the optim optimiser.}

\item{lower_initialisation}{Lower bound for the uniform random variable used
to generate initial values.}

\item{upper_initialisation}{Upper bound for the uniform random variable used
to generate initial values.}

\item{n_restarts}{Number of different initialisation values to test.}
}
\value{
A \eqn{4 x n_samples} matrix. Hyperparameters are ordered in the rows
  as (phylogenetic noise, phylogenetic length-scale, non-phylogenetic noise).
  Each column is a separate sample.
}
\description{
Type II Maximum Likelihood Estimation of hyperparameters for n_samples
realisations of a Phylogenetic Ornstein-Uhlenbeck Process with given
hyperparameters
}
