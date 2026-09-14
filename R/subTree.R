#' @include ConvenientTblTree.R
NULL

#' Extract a subtree
#'
#' Takes a `tbl_tree` object representing a tree, and extracts a subtree at a
#' given node.
#'
#' @param tree_tibble A [`ConvenientTblTree`] object.
#' @param node The ID number of the node to base the subtree, or a `FocalClade`
#'        object.
#' @param ... Other arguments (currently ignored)
#'
#' @importFrom tidytree as_tibble as.treedata tree_subset
#'
#' @return A subtree represented in the `tbl_tree` class of the _tidytree_
#' package.  This tibble has an extra `node.orig` column providing the old node
#' IDs, because the node IDs in the object are renamed during the extraction
#' of the subtree.
#'
#' @examples
#' (subtree <- subTree(Halo_Tree, 9))
#' visualizeTree(subtree, subtree$node.orig)
#'
#' @author Charles Plessy
#' @family Focal clade functions
#'
#' @export
subTree <- new_generic("subTree", c("tree_tibble", "node"))

method(subTree, list(ConvenientTblTree, class_numeric)) <- function(tree_tibble, node) {
  tree_tibble$node.orig <- tree_tibble$node
  tree_tibble |> as.treedata() |> tree_subset(node, levels_back = 0) |> as_tibble() |> suppressMessages() |> makeItConvenient()
}
