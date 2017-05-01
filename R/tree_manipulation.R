#' Tree expansion
#'
#' Expand a phylogenetic tree by adding extra nodes to each existing tip of the tree.
#'
#' @inheritParams pou_covariance
#' @param n_nodes The number of nodes to be added to each tip of the tree.
#' @param evolutionary_time The distance from the root of the nodes being added to each of the n_nodes tips.
#'
#' @return An object of class 'phylo' with a similar structure to the input \eqn{phylogenetic_tree} but with \eqn{n_nodes} node at each tip.
expand_tree <- function(phylogenetic_tree, n_nodes, evolutionary_time = 0){
  expanded <- phylogenetic_tree
  tips_for_expansion <- expanded$tip.label

  for(i in tips_for_expansion){
    expansion <- paste("(",paste(rep(paste0(i, ':', evolutionary_time), n_nodes), collapse = ", ")  ,");" )
    expansion_tree <- read.tree(text = expansion)

    place <- which(expanded$tip.label == i)

    expanded <-  bind.tree(expanded, expansion_tree, place)
  }
  return(expanded)
}
