#' @include FocalClade.R
NULL

CladeBars <- new_class("CladeBars", properties = list(
  clades = FocalCladeList,
  offset = class_numeric,
  width  = class_numeric
))

#' Focal clades as bars
#'
#' Mark focal clades with colored vertical bars to the right of the tips of a
#' tree plot, instead of highlighting them with colored boxes behind the
#' branches.  This avoids color clashes with plots whose branches are colored,
#' such as the ones of [visualizeKaryotype()].
#'
#' The bars are drawn with the `fill` aesthetic, so their legend is separate
#' from the `colour` legend of the branches.  Nested clades are drawn in
#' separate columns, the largest one closest to the tree.
#'
#' @param clades A [`FocalCladeList`] or a single [`FocalClade`] object.
#' @param offset Space between the tips and the bars, as a fraction of the
#'        tree's width.
#' @param width Width of the bars, as a fraction of the tree's width.
#'
#' @returns An object to add with `+` to a tree plot made with
#' [visualizeTree()] or [visualizeKaryotype()].  Clades whose genomes are not
#' in the plotted tree are skipped.
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Focal clade functions
#' @family Plotting functions
#'
#' @examples
#' visualizeTree(Halo_Tree) + cladeBars(Halo_FocalClades)
#'
#' @importFrom ggplot2 aes geom_rect scale_fill_manual
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export

cladeBars <- function(clades, offset = 0.1, width = 0.02) {
  if (S7_inherits(clades, FocalClade)) clades <- FocalCladeList(clades)
  CladeBars(clades = clades, offset = offset, width = width)
}

method(`+`, list(PlottedTree, CladeBars)) <- function(e1, e2) {
  tips  <- e1$data[e1$data$isTip, ]
  xmax  <- max(e1$data$x)
  bars  <- lapply(e2@clades, \(clade) {
    y <- tips$y[tips$label %in% clade@genomeIDs]
    if (length(y) == 0) return(NULL)
    data.frame(clade = clade@displayName, color = clade@color,
               ymin = min(y) - 0.4, ymax = max(y) + 0.4)
  }) |> do.call(what = rbind)
  bars$clade <- factor(bars$clade, levels = bars$clade)
  bars$xmin  <- xmax * (1 + e2@offset + 1.5 * e2@width * barColumns(bars$ymin, bars$ymax))
  bars$xmax  <- bars$xmin + xmax * e2@width
  e1 +
    geom_rect(data = bars,
              aes(xmin = .data$xmin, xmax = .data$xmax,
                  ymin = .data$ymin, ymax = .data$ymax, fill = .data$clade),
              inherit.aes = FALSE) +
    scale_fill_manual(values = setNames(bars$color, bars$clade), name = "Clade")
}

# Assign overlapping bars (nested clades) to successive columns, starting from
# the largest clade, so that each bar takes the first column where it does not
# overlap another one.  Returns the column index of each bar, starting at 0.
barColumns <- function(ymin, ymax) {
  column <- rep(NA_integer_, length(ymin))
  for (i in order(ymax - ymin, decreasing = TRUE)) {
    taken <- column[!is.na(column) & ymin <= ymax[i] & ymax >= ymin[i]]
    column[i] <- min(setdiff(0:length(ymin), taken))
  }
  column
}
