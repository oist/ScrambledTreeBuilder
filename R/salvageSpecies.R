#' Species to salvage when subsampling
#'
#' Find the species whose position in the tree is the most at risk of being
#' wrong, so that they can be kept by [subsampleSpecies()]: in a subsample,
#' a misplaced lineage that is missing can not be detected.
#'
#' Two kinds of species are salvaged:
#'
#' - The species on the `long` longest terminal branches.  Their closest
#'   relatives in the data set are far, so their placement rests on few
#'   shared changes and is prone to errors such as long-branch attraction.
#' - For each internal node deeper than `depth` whose bootstrap support is
#'   below `support`, at least one species on each side of the node, so that
#'   every lineage whose position is uncertain stays represented whatever
#'   the true topology is.  Nodes are visited from the deepest one, and on
#'   each side the species with the longest terminal branch is chosen,
#'   unless a species of that side is already salvaged or in `keep`.
#'
#' @param tree A [`stbTree`] object, with a `branch.length` column and, to
#'        salvage the sides of weakly supported nodes, a `bootstrap` column
#'        (see the `n_bootstrap` option of [makeTidyTree()]).
#' @param value Name of a column of the tree containing the values of the
#'        internal nodes, such as a pairwise distance, or a vector of these
#'        values.  Used to measure the depth of nodes.
#' @param long Number of species to salvage for their long terminal
#'        branches.
#' @param support Bootstrap support below which a node is weakly supported,
#'        in the units of the `bootstrap` column.  `NULL` to salvage only
#'        long branches.
#' @param depth Node value above which weakly supported nodes are
#'        considered.  By default, half of the largest node value: shallow
#'        nodes, between closely related genomes, are often weakly supported
#'        without consequence.
#' @param keep Species already kept, which count as representatives of the
#'        sides of weakly supported nodes.
#'
#' @returns A data frame with one row per salvaged species, with columns
#' `species` and `reason`.  Pass `c(keep, salvaged$species)` to the `keep`
#' argument of [subsampleSpecies()].
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Functions for trees
#' @seealso [subsampleSpecies()]
#'
#' @examples
#' salvageSpecies(Halo_Tree, "Percent_difference", long = 2, support = NULL)
#'
#' @export

salvageSpecies <- function(tree, value, long = 5, support = 70, depth = NULL, keep = NULL) {
  valueName <- if (length(value) == 1 && is.character(value)) value else "value"
  if (length(value) == 1 && is.character(value)) {
    if (is.null(tree[[value]])) stop("Column ", dQuote(value), " not found in the tree.")
    value <- tree[[value]]
  }
  if (!"branch.length" %in% names(tree)) stop("The tree has no branch.length column.")
  tips   <- tree$isTip
  tipBl  <- setNames(tree$branch.length[tips], tree$label[tips])
  out    <- data.frame(species = character(), reason = character())

  if (long > 0) {
    longest <- head(names(sort(tipBl, decreasing = TRUE)), long)
    out <- rbind(out, data.frame(species = longest,
                                 reason  = sprintf("long terminal branch (%.3g)", tipBl[longest])))
  }

  if (!is.null(support)) {
    if (!"bootstrap" %in% names(tree) || all(is.na(tree$bootstrap[!tips]))) {
      warning("No bootstrap support in the tree: weakly supported nodes are not salvaged.")
    } else {
      inner <- tree$node[!tips]
      v <- value[match(inner, tree$node)]
      b <- tree$bootstrap[match(inner, tree$node)]
      if (is.null(depth)) depth <- max(v, na.rm = TRUE) / 2
      weak <- which(!is.na(b) & b < support & v >= depth)
      chosen <- c(keep, out$species)
      for (i in weak[order(-v[weak])]) {
        for (side in childSpecies(tree, inner[i])) {
          if (any(side %in% chosen)) next
          pick <- side[which.max(tipBl[side])]
          chosen <- c(chosen, pick)
          out <- rbind(out, data.frame(species = pick,
                                       reason  = sprintf("weakly supported node %d (support %g, %s %.3g)",
                                                         inner[i], b[i], valueName, v[i])))
        }
      }
    }
  }
  rownames(out) <- NULL
  out
}
