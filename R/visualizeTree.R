#' Plot a phylogenetic tree and its node data
#'
#' @param tree A [`stbTree`] object.
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

visualizeTree <- function(tree, value="node", valueround = 2, outerlabelsize = 0.25, innerlabelsize = 3, ynudge = 0, xnudge = 0) {
  noLegend <- FALSE
  if (length(value) == 1 && value == "node") noLegend <- TRUE
  if (length(value) == 1 && is.character(value))
    value <- tree[ , value, drop = TRUE]
  # Build step by step for better use of suppressMessages
  gg <- ggtree::ggtree(tidytree::as.treedata(tree))
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

# Label the nodes of a plotted tree with values, or mark them with points
# colored by the values.  Shared by visualizeTree() and visualizeKaryotype().
#' @importFrom ggplot2 aes geom_label geom_point position_nudge unit
#' @importFrom rlang .data
#' @noRd
addValuesToTree <- function (ggtree, value, valueround = 2, outerlabelsize = 0.25, innerlabelsize = 3, ynudge = 0, xnudge = 0, points = FALSE) {
  if (isTRUE(points)) {
    # Only nodes with a value: unlike labels, points with NA values are drawn.
    withValues <- function(d) {
      d$value <- value
      d[!is.na(d$value), ]
    }
    return(ggtree + geom_point( data = withValues
                              , aes(color = .data$value)
                              , size = innerlabelsize
                              , position = position_nudge(x = xnudge, y = ynudge)))
  }
  ggtree + geom_label( aes( label = round(value, digits = valueround)
                          , color = value)
                     , linewidth = outerlabelsize
                     , linetype = "solid" # Not inherited from dashed branches.
                     , size = innerlabelsize
                     , na.rm = TRUE
                     , label.padding = unit(0.15, "lines")
                     , nudge_y = ynudge, nudge_x = xnudge)
}
