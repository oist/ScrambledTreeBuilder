#' Plot a phylogenetic tree and its node data
#'
#' Plot a tree, label its internal nodes with values, and optionally color its
#' branches by a species trait such as the number of chromosomes.
#'
#' @section Species traits:
#' A trait is a value that belongs to each species, not to pairs of species,
#' for instance a chromosome number.  Store it in a column of the tree, with
#' values on the tips; values on the internal nodes are ignored.  The values
#' of the internal nodes are reconstructed with [ancestralStates()], by
#' parsimony unless another method is chosen with the `...` options.  Each
#' branch takes the color of the node it leads to, and the tip values are
#' printed next to the species names.  Branches leading to nodes whose state
#' is ambiguous are dashed and light grey.  Values are treated as discrete:
#' the most common one is dark grey and each of the others has its own color
#' from the colorblind-friendly Okabe-Ito palette (or from a larger palette
#' when there are more than seven).  Node values then have their own color
#' scale.  To mark focal clades without color clashes, add them as bars with
#' [cladeBars()].
#'
#' @param tree A [`stbTree`] object.
#' @param value Values to label on internal nodes of the tree, or name of a
#'        column in the tree object.  By default, the node IDs, without
#'        legend.  `NULL` labels nothing.
#' @param valueround Number of decimals to round values.
#' @param outerlabelsize Size of label border.
#' @param innerlabelsize Overall size of label (or of points).
#' @param xnudge,ynudge Adjust horizontal or vertical position of labels
#'        (useful when plotting multiple labels).
#' @param points Mark the nodes with points colored by value instead of
#'        labels, for large trees.
#' @param trait Name of a column of the tree containing a species trait on
#'        its tips, to color the branches (see below).  By default, no trait
#'        is shown.
#' @param offset Space between the tips and the column of trait values, as a
#'        fraction of the tree's width.
#' @param background The trait value drawn in grey, so that the other values
#'        stand out.  By default, the most common one among the tips.
#' @param colors A vector of colors named after the trait values, to replace
#'        the default palette.
#' @param axis Draw an axis of pairwise distances: for each node, the
#'        distance between the species on each side of it, assuming the tree
#'        was built by [makeTidyTree()] (UPGMA), so that the height of a node
#'        is half that distance.  `"number"` (default) labels the distances
#'        as they are, `"percent"` multiplies them by 100 and adds a percent
#'        sign, for distances that are fractions, and `"none"` draws no axis.
#'        Unambiguous abbreviations are accepted.  The axis has no title: add
#'        one with [ggplot2::labs()].
#' @param ... Options passed to [ancestralStates()] to choose the method of
#'        reconstruction of the `trait` (`method`, `model`, `threshold`).
#'
#' @return A `ggtree` plot, to which focal clades can be added (see
#' [focalClade()] and [cladeBars()]).
#'
#' @export
#'
#' @author Noa Brenner
#' @author Charles Plessy
#'
#' @family Focal clade functions
#' @family Plotting functions
#' @seealso [ancestralStates()]
#'
#' @examples
#' visualizeTree(Halo_Tree)
#' visualizeTree(Halo_Tree, value = Halo_Tree$Scrambling_index)
#' visualizeTree(Halo_Tree, "Scrambling_index") +  # same
#'   ggplot2::labs(title = "Scrambling index", x = "Pairwise percent difference")
#'
#' # Pretend chromosome numbers (the real ones are all 1), stored on the tips.
#' chr <- c(Halobacterium_litoreum = 3, Halobacterium_noricense = 3,
#'          Halobacterium_salinarum = 3, Salarchaeum_japonicum  = 5,
#'          Haloferax_mediterranei  = 2, Haloferax_volcanii     = 2)
#' tree <- Halo_Tree
#' tree$ChromNumber <- chr[tree$label]
#' visualizeTree(tree, value = NULL, trait = "ChromNumber") + cladeBars(Halo_FocalClades)
#'
#' # Node values have their own colour scale, which can be replaced.
#' visualizeTree(tree, "Scrambling_index", trait = "ChromNumber", points = TRUE) +
#'   ggplot2::scale_colour_viridis_c(name = "Scrambling index", option = "magma")
#'
#' @importFrom ggplot2 guides labs
#' @importFrom ggnewscale new_scale_colour
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom tidytree as.treedata

visualizeTree <- function(tree, value = "node", valueround = 2, outerlabelsize = 0.25,
                          innerlabelsize = 3, ynudge = 0, xnudge = 0, points = FALSE,
                          trait = NULL, offset = 0.05, background = NULL, colors = NULL,
                          axis = c("number", "percent", "none"), ...) {
  axis <- match.arg(axis)
  nodeIDs   <- identical(value, "node")
  valueName <- if (length(value) == 1 && is.character(value)) value else "value"
  if (length(value) == 1 && is.character(value))
    value <- tree[ , value, drop = TRUE]
  if (is.null(trait)) {
    gg <- ggtree(as.treedata(tree))
    # Build step by step for better use of suppressMessages
    suppressMessages(
      gg <- gg + geom_tiplab(as_ylab = TRUE)
    ) # Scale for y is already present.
  } else {
    gg <- traitTree(tree, trait, offset, background, colors, ...)
    if (!is.null(value)) gg <- gg + new_scale_colour()
  }
  if (!is.null(value)) {
    gg <- addValuesToTree(gg, value = value, valueround = valueround,
                          outerlabelsize = outerlabelsize, innerlabelsize = innerlabelsize,
                          ynudge = ynudge, xnudge = xnudge, points = points) +
      labs(colour = valueName)
    if (nodeIDs) gg <- gg + guides(colour = "none")
  }
  if (axis != "none") gg <- addPairwiseAxis(gg, percent = axis == "percent")
  gg
}

# Plot a tree with branches colored by a species trait reconstructed from the
# values on its tips, and the tip values in a column next to the tip labels.
#' @importFrom ggplot2 aes geom_text guide_legend guides labs
#' @importFrom ggplot2 scale_colour_manual scale_linetype_manual
#' @importFrom grDevices hcl.colors
#' @importFrom stats na.omit setNames
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom rlang .data
#' @importFrom tidytree as.treedata
#' @noRd
traitTree <- function(tree, trait, offset = 0.05, background = NULL, colors = NULL, ...) {
  if (!is.character(trait) || length(trait) != 1)
    stop("trait must be the name of a column of the tree.")
  if (is.null(tree[[trait]]))
    stop("Column ", dQuote(trait), " not found in the tree.")
  tips   <- tree$isTip
  taxons <- setNames(data.frame(tree[[trait]][tips], row.names = tree$label[tips]), trait)
  tree   <- ancestralStates(tree, taxons, trait, ...)
  states <- sort(unique(na.omit(tree[[trait]])))
  if (is.null(colors)) {
    if (is.null(background)) {
      tipStates  <- tree[[trait]][tree$isTip]
      background <- states[which.max(tabulate(match(tipStates, states)))]
    }
    others <- setdiff(states, background)
    palette <- if (length(others) <= length(traitPalette)) traitPalette
               else hcl.colors(length(others), "Dark 3")
    colors <- setNames(c("grey30", palette[seq_along(others)]),
                       as.character(c(background, others)))
  }
  # Legend keys: one solid line per state, then a dashed one if ambiguous.
  keyTypes <- c(rep("solid", length(states)), if (anyNA(tree[[trait]])) "dashed")
  valueColumn <- function(d) {
    d <- d[d$isTip, ]
    d$x <- max(d$x) * (1 + offset)
    d
  }
  gg <- ggtree(as.treedata(tree), aes(color    = factor(.data[[trait]]),
                                      linetype = is.na(.data[[trait]])))
  # Build step by step for better use of suppressMessages
  suppressMessages(
    gg <- gg + geom_tiplab(as_ylab = TRUE)
  ) # Scale for y is already present.
  gg +
    geom_text(data = valueColumn,
              aes(x = .data$x, y = .data$y, label = .data[[trait]]),
              inherit.aes = FALSE, hjust = 0, size = 3) +
    scale_colour_manual(values = colors, na.value = "grey70",
                        labels = \(x) ifelse(is.na(x), "ambiguous", x)) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    guides(colour = guide_legend(override.aes = list(linetype = keyTypes))) +
    labs(colour = trait)
}

# Okabe-Ito colors, without black and grey (used for the background and
# ambiguous states), and with the low-contrast yellow last.
traitPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")

# Add an axis of pairwise distances, increasing from the tips (0) towards the
# root.  With UPGMA, the height of a node is half the distance between the
# species on each side of it.  With percent = TRUE, the distances are
# fractions to be labelled as percentages.
#' @importFrom ggplot2 element_line element_text scale_x_continuous theme
#' @noRd
addPairwiseAxis <- function(gg, percent = FALSE) {
  xmax   <- max(gg$data$x[gg$data$isTip])
  breaks <- pretty(c(0, 2 * xmax))
  breaks <- breaks[breaks <= 2 * xmax]
  labels <- if (percent) paste0(100 * breaks, "%") else breaks
  suppressMessages( # Scale for x is already present.
    gg + scale_x_continuous(breaks = xmax - breaks / 2, labels = labels)
  ) + theme(axis.line.x  = element_line(),
            axis.ticks.x = element_line(),
            axis.text.x  = element_text())
}

# Label the nodes of a plotted tree with values, or mark them with points
# colored by the values.
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
