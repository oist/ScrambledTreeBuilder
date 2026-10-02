#' Reconstruct ancestral states of a species trait
#'
#' Given a value for each species (tip) of a tree, such as its chromosome
#' number, infer the value of each internal node with Fitch parsimony, that
#' is, with the smallest number of changes along the branches.
#'
#' The values are treated as discrete states: 12 is as different from 13 as
#' from 24.  The first pass of the algorithm collects, from the tips to the
#' root, the set of states each node could have.  The second pass, from the
#' root to the tips, narrows each node's set to the states it shares with its
#' parent, if any.  Nodes left with more than one possible state are
#' ambiguous and get `NA`.  Tips keep their original values.
#'
#' @param tree A [`stbTree`] object.
#' @param values Either a named vector of values whose names are the tip labels
#'        of the tree, or a taxon table (a data frame with one row per species
#'        and the tip labels as row names).
#' @param column When `values` is a data frame, the name of the column
#'        containing the values.
#'
#' @returns The `tree` with a new column, named after `column` or `"state"`,
#' containing the tip values and the reconstructed internal node values, with
#' `NA` for ambiguous nodes.  Missing tip values are treated as unknown (any
#' state) during reconstruction, with a warning.
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Functions for trees
#' @seealso [visualizeKaryotype()]
#'
#' @examples
#' # Pretend chromosome numbers (the real ones are all 1).
#' chr <- c(Halobacterium_litoreum = 3, Halobacterium_noricense = 3,
#'          Halobacterium_salinarum = 3, Salarchaeum_japonicum  = 5,
#'          Haloferax_mediterranei  = 2, Haloferax_volcanii     = 2)
#' ancestralStates(Halo_Tree, chr)
#'
#' @importFrom stats na.omit setNames
#' @export

ancestralStates <- function(tree, values, column = NULL) {
  colname <- if (is.null(column)) "state" else column
  if (is.data.frame(values)) {
    if (is.null(column) || is.null(values[[column]]))
      stop("Column ", dQuote(column), " not found in the taxon table.")
    values <- setNames(values[[column]], rownames(values))
  }
  if (is.null(names(values)))
    stop("Values must be named after the tip labels of the tree.")

  tips <- tree$node[tree$isTip]
  tipValues <- unname(values[tree$label[tips]])
  missing <- tree$label[tips][is.na(tipValues)]
  if (length(missing) > 0)
    warning("No value for: ", paste(missing, collapse = ", "), ".  Treating as unknown.")

  states <- sort(unique(na.omit(tipValues)))
  root   <- tree$node[tree$parent == tree$node]

  # Possible states of each node, as a logical matrix (nodes x states).
  sets <- matrix(FALSE, nrow = max(tree$node), ncol = length(states))
  sets[cbind(tips, match(tipValues, states))[!is.na(tipValues), , drop = FALSE]] <- TRUE
  sets[tips[is.na(tipValues)], ] <- TRUE

  children <- function(node) tree$node[tree$parent == node & tree$node != node]

  # First pass, from the tips to the root: keep the states shared by the most
  # children (their intersection for a binary node, or else their union).
  upPass <- function(node) {
    if (node %in% tips) return(invisible())
    kids <- children(node)
    lapply(kids, upPass)
    counts <- colSums(sets[kids, , drop = FALSE])
    sets[node, ] <<- counts == max(counts)
  }
  upPass(root)

  # Second pass, from the root to the tips: narrow each node's states to the
  # ones it shares with its parent.
  downPass <- function(node) {
    for (kid in setdiff(children(node), tips)) {
      shared <- sets[kid, ] & sets[node, ]
      if (any(shared)) sets[kid, ] <<- shared
      downPass(kid)
    }
  }
  downPass(root)

  resolved <- rowSums(sets) == 1
  result <- states[NA_integer_]
  result[resolved[tree$node]] <- states[apply(sets[tree$node[resolved[tree$node]], , drop = FALSE], 1, which)]
  result[tree$isTip] <- tipValues[match(tree$node[tree$isTip], tips)]
  tree[[colname]] <- result
  tree
}
