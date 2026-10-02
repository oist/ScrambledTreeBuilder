#' Plot karyotype evolution on a tree
#'
#' Color the branches of a tree by chromosome number, so that lineages keep
#' the same color until their karyotype changes, and print the chromosome
#' numbers next to the species names.
#'
#' The chromosome numbers of internal nodes are reconstructed with
#' [ancestralStates()], by parsimony unless another method is chosen with the
#' `...` options.  Each branch takes the color of the node it leads to.
#' Branches leading to nodes whose state is ambiguous are dashed and light
#' grey.  Chromosome numbers are treated as discrete values: the most common
#' one is dark grey and each of the others has its own color from the
#' colorblind-friendly Okabe-Ito palette (or from a larger palette when there
#' are more than seven).
#'
#' To mark focal clades without color clashes, add them as bars with
#' [cladeBars()].
#'
#' @param tree A [`stbTree`] object.
#' @param taxons A taxon table: a data frame with one row per species and the
#'        tip labels of the tree as row names.
#' @param column Name of the column of `taxons` containing the chromosome
#'        numbers.
#' @param offset Space between the tips and the column of numbers, as a
#'        fraction of the tree's width.
#' @param background The chromosome number drawn in grey, so that the other
#'        numbers stand out.  By default, the most common one among the tips.
#' @param colors A vector of colors named after the chromosome numbers, to
#'        replace the default palette.
#' @param value Values to label on the internal nodes of the tree, or name of
#'        a column in the tree object, like in [visualizeTree()].  By default,
#'        no values are shown.
#' @param valueround Number of decimals to round values.
#' @param points Mark the nodes with points colored by value instead of
#'        labels, for large trees.
#' @param axis Draw an axis of pairwise distances: for each node, the
#'        distance between the species on each side of it, assuming the tree
#'        was built by [makeTidyTree()] (UPGMA), so that the height of a node
#'        is half that distance.  The axis has no title: add one with
#'        [ggplot2::xlab()].
#' @param ... Options passed to [ancestralStates()] to choose the method of
#'        reconstruction (`method`, `model`, `threshold`).
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
#' visualizeKaryotype(Halo_Tree, taxons) + cladeBars(Halo_FocalClades)
#'
#' # Overlay node values, with an axis of the distances used to build the tree.
#' visualizeKaryotype(Halo_Tree, taxons, value = "Scrambling_index", axis = TRUE) +
#'   ggplot2::xlab("Pairwise percent difference")
#'
#' # Node values have their own colour scale, which can be replaced.
#' visualizeKaryotype(Halo_Tree, taxons, value = "Scrambling_index", points = TRUE) +
#'   ggplot2::scale_colour_viridis_c(name = "Scrambling index", option = "magma")
#'
#' @importFrom ggplot2 aes geom_text guide_legend guides labs
#' @importFrom ggplot2 scale_colour_continuous scale_colour_manual scale_linetype_manual
#' @importFrom grDevices hcl.colors
#' @importFrom stats na.omit setNames
#' @importFrom ggnewscale new_scale_colour
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom rlang .data
#' @importFrom tidytree as.treedata
#' @export

visualizeKaryotype <- function(tree, taxons, column = "ChromNumber", offset = 0.05,
                               background = NULL, colors = NULL, value = NULL,
                               valueround = 2, points = FALSE, axis = FALSE, ...) {
  tree <- ancestralStates(tree, taxons, column, ...)
  valueName <- if (length(value) == 1 && is.character(value)) value else "value"
  if (length(value) == 1 && is.character(value))
    value <- tree[ , value, drop = TRUE]
  states <- sort(unique(na.omit(tree[[column]])))
  if (is.null(colors)) {
    if (is.null(background)) {
      tipStates  <- tree[[column]][tree$isTip]
      background <- states[which.max(tabulate(match(tipStates, states)))]
    }
    others <- setdiff(states, background)
    palette <- if (length(others) <= length(karyotypePalette)) karyotypePalette
               else hcl.colors(length(others), "Dark 3")
    colors <- setNames(c("grey30", palette[seq_along(others)]),
                       as.character(c(background, others)))
  }
  # Legend keys: one solid line per state, then a dashed one if ambiguous.
  keyTypes <- c(rep("solid", length(states)), if (anyNA(tree[[column]])) "dashed")
  numberColumn <- function(d) {
    d <- d[d$isTip, ]
    d$x <- max(d$x) * (1 + offset)
    d
  }
  gg <- ggtree(as.treedata(tree), aes(color    = factor(.data[[column]]),
                                      linetype = is.na(.data[[column]])))
  # Build step by step for better use of suppressMessages
  suppressMessages(
    gg <- gg + geom_tiplab(as_ylab = TRUE)
  ) # Scale for y is already present.
  gg <- gg +
    geom_text(data = numberColumn,
              aes(x = .data$x, y = .data$y, label = .data[[column]]),
              inherit.aes = FALSE, hjust = 0, size = 3) +
    scale_colour_manual(values = colors, na.value = "grey70",
                        labels = \(x) ifelse(is.na(x), "ambiguous", x)) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    guides(colour = guide_legend(override.aes = list(linetype = keyTypes))) +
    labs(colour = column)
  if (!is.null(value)) {
    gg <- addValuesToTree(gg + new_scale_colour(), value = value,
                          valueround = valueround, points = points) +
      scale_colour_continuous(name = valueName)
  }
  if (isTRUE(axis)) gg <- addPairwiseAxis(gg)
  gg
}

# Add an axis of pairwise distances, increasing from the tips (0) towards the
# root.  With UPGMA, the height of a node is half the distance between the
# species on each side of it.
#' @importFrom ggplot2 element_line element_text scale_x_continuous theme
#' @noRd
addPairwiseAxis <- function(gg) {
  xmax   <- max(gg$data$x[gg$data$isTip])
  breaks <- pretty(c(0, 2 * xmax))
  breaks <- breaks[breaks <= 2 * xmax]
  suppressMessages( # Scale for x is already present.
    gg + scale_x_continuous(breaks = xmax - breaks / 2, labels = breaks)
  ) + theme(axis.line.x  = element_line(),
            axis.ticks.x = element_line(),
            axis.text.x  = element_text())
}

# Okabe-Ito colors, without black and grey (used for the background and
# ambiguous states), and with the low-contrast yellow last.
karyotypePalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
