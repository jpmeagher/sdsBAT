#' Tree expansion
#'
#' Expand a phylogenetic tree by adding extra nodes to each existing tip of the tree.
#'
#' @inheritParams pou_covariance
#' @param n_nodes The number of nodes to be added to each tip of the tree.
#' @param evolutionary_time The distance from the root of the nodes being added to each of the n_nodes tips.
#'
#' @return An object of class 'phylo' with a similar structure to the input \eqn{phylogenetic_tree} but with \eqn{n_nodes} node at each tip.
#'
#' @export
expand_tree <- function(phylogenetic_tree, n_nodes, evolutionary_time = 0){
  if(attributes(phylogenetic_tree)$class != 'phylo'){
    stop("phylogenetic tree must be an object of class 'phylo'.")
  }
  if(length(n_nodes) != 1){
    stop("Input number of nodes to be added to each tip")
  }
  if(length(evolutionary_time) != 1){
    stop("Input the distance from existing tips to added nodes")
  }

  expanded <- phylogenetic_tree
  tips_for_expansion <- expanded$tip.label

  for(i in tips_for_expansion){
    expansion <- paste("(",paste(rep(paste0(i, ':', evolutionary_time), n_nodes), collapse = ", ")  ,");" )
    expansion_tree <- ape::read.tree(text = expansion)

    place <- which(expanded$tip.label == i)

    expanded <-  ape::bind.tree(expanded, expansion_tree, place)
  }
  return(expanded)
}
