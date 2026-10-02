#' Plot karyotype evolution on a tree
#'
#' Color the branches of a tree by chromosome number, so that lineages keep
#' the same color until their karyotype changes, and print the chromosome
#' numbers next to the species names.
#'
#' The chromosome numbers of internal nodes are reconstructed with
#' [ancestralStates()].  Each branch takes the color of the node it leads to.
#' Branches leading to nodes whose state is ambiguous are grey.  Chromosome
#' numbers are treated as discrete values, each with its own color.
#'
#' @param tree A [`stbTree`] object.
#' @param taxons A taxon table: a data frame with one row per species and the
#'        tip labels of the tree as row names.
#' @param column Name of the column of `taxons` containing the chromosome
#'        numbers.
#' @param offset Space between the tips and the column of numbers, as a
#'        fraction of the tree's width.
#'
#' @returns A `ggtree` plot, to which focal clades can be added (see
#' [focalClade()]).
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Plotting functions
#' @seealso [ancestralStates()], [visualizeTree()]
#'
#' @examples
#' # Pretend chromosome numbers (the real ones are all 1).
#' taxons <- data.frame(row.names = Halo_Tree$label[Halo_Tree$isTip],
#'                      ChromNumber = c(3, 3, 3, 2, 2, 5))
#' visualizeKaryotype(Halo_Tree, taxons)
#' visualizeKaryotype(Halo_Tree, taxons) + Halo_FocalClades
#'
#' @importFrom ggplot2 aes geom_text labs scale_colour_discrete
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom rlang .data
#' @importFrom tidytree as.treedata
#' @export

visualizeKaryotype <- function(tree, taxons, column = "ChromNumber", offset = 0.05) {
  tree <- ancestralStates(tree, taxons, column)
  numberColumn <- function(d) {
    d <- d[d$isTip, ]
    d$x <- max(d$x) * (1 + offset)
    d
  }
  gg <- ggtree(as.treedata(tree), aes(color = factor(.data[[column]])))
  # Build step by step for better use of suppressMessages
  suppressMessages(
    gg <- gg + geom_tiplab(as_ylab = TRUE)
  ) # Scale for y is already present.
  gg +
    geom_text(data = numberColumn,
              aes(x = .data$x, y = .data$y, label = .data[[column]]),
              inherit.aes = FALSE, hjust = 0, size = 3) +
    scale_colour_discrete(labels = \(x) ifelse(is.na(x), "ambiguous", x)) +
    labs(colour = column)
}
