#' @include stbTree.R
NULL

#' Make a tidy tree
#'
#' Given a distance matrix, computes a cladogram with the [`stats::hclust`]
#' function, and returns it as a [`tibble::tibble`] formatted for the `tidytree`
#' package.
#'
#' @param pairwise_matrix a distance matrix
#' @param n_bootstrap Integer number of bootstrap replicates (default 0, no bootstrap).
#'
#' @return A `tbl_tree` object.  Bootstrap values will be `NA` if no boostrap is
#' done.
#'
#' @importFrom stats hclust as.dist
#'
#' @examples
#' makeTidyTree((Halo_PercentDiff/2 + t(Halo_PercentDiff)/2))
#'
#' @export

makeTidyTree <- function(pairwise_matrix, n_bootstrap = 0) {
  if (! all(pairwise_matrix == t(pairwise_matrix)))
    warning("Input matrix is not symmetric!")
  if (n_bootstrap == 0) {
    HClust <- hclust(as.dist(pairwise_matrix), method = "average")
    tree <- tidytree::as_tibble(tidytree::as.phylo(HClust))
  } else {
    bt <- residualBootstrapTree(pairwise_matrix, n_bootstrap = n_bootstrap)
    tree <- tidytree::as_tibble(tidytree::as.phylo(bt$tree))
  }
  makeItConvenient(tree)
}

#' Add needed columns to be convenient
#'
#' The `stbTree` S7 class is only a draft and lacks validation
#' functions that ensure that an object provides everything expected from
#' such a tree, in particular the `isTip` and `y` columns that can become
#' invalid upon subsetting.
#'
#' This private function (re)computes `isTip` and `y`.
#'
#' @param tree A [`stbTree`] object.
#'
#' @returns A [`stbTree`] object with valid `isTip` and `y`
#' information.

makeItConvenient <- function(tree) {
  plotted_tree <- visualizeTree(tree)$data # where tree gains isTip, x, y, branch and angle columns
  tree$isTip <- plotted_tree$isTip
  tree$y     <- plotted_tree$y
  tree$bootstrap <- tree$label |> as.numeric() |> suppressWarnings() # I know I discard character values
  tree$label[!tree$isTip] <- NA
  tree
}
