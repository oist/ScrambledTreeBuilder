#' Plot a phylogenetic tree and its node data
#'
#' @param your_tibble A [`stbTree`] object.
#' @param value Tibble value to label on internal nodes of the tree, or name
#'        of a column in the tree object.
#' @param valueround Number of integers to round value.
#' @param outerlabelsize Size of label border.
#' @param innerlabelsize Overall size of label.
#' @param xnudge,ynudge Adjust horizontal or vertical position of labels
#'        (useful when plotting multiple labels).
#'
#' @return A phylogenetic tree showcasing phylogeny of species, additional node
#' values, or just node IDs.  No legend is show when displaying node IDs.
#'
#' @export
#'
#' @author Noa Brenner
#' @author Charles Plessy
#'
#' @family Focal clade functions
#' @family Plotting functions
#'
#' @examples
#' visualizeTree(Halo_Tree)
#' visualizeTree(Halo_Tree, value = Halo_Tree$Scrambling_index)
#' visualizeTree(Halo_Tree, "Scrambling_index") # same

visualizeTree <- function(your_tibble, value="node", valueround = 2, outerlabelsize = 0.25, innerlabelsize = 3, ynudge = 0, xnudge = 0) {
  addValuesToTree <- function (ggtree, value, valueround = 2, outerlabelsize = 0.25, innerlabelsize = 3, ynudge = 0, xnudge = 0) {
    ggtree + ggtree::geom_label( ggtree::aes( label = round(value, digits = valueround)
                                            , color = value)
                               , label.size = outerlabelsize
                               , size = innerlabelsize
                               , na.rm = TRUE
                               , label.padding = ggtree::unit(0.15, "lines")
                               , nudge_y = ynudge, nudge_x = xnudge)
  }
  noLegend <- FALSE
  if (length(value) == 1 && value == "node") noLegend <- TRUE
  if (length(value) == 1 && is.character(value))
    value <- your_tibble[ , value, drop = TRUE]
  # Build step by step for better use of suppressMessages
  gg <- ggtree::ggtree(tidytree::as.treedata(your_tibble))
  suppressMessages(
    gg <- gg + ggtree::geom_tiplab(as_ylab=TRUE)
  ) # Scale for y is already present.
  gg <- addValuesToTree(gg, value = value, valueround = valueround, outerlabelsize = outerlabelsize, innerlabelsize = innerlabelsize, ynudge = ynudge, xnudge = xnudge)
  if (noLegend) {
    gg + ggplot2::theme(legend.position = "none")
  } else {
    gg
  }
}
