#' Reconstruct ancestral states of a species trait
#'
#' Given a value for each species (tip) of a tree, such as its chromosome
#' number, infer the value of each internal node, either with Fitch parsimony
#' (the smallest number of changes along the branches) or by maximum
#' likelihood with [ape::ace()].
#'
#' @section Parsimony:
#' The values are treated as discrete states: 12 is as different from 13 as
#' from 24.  The first pass of the algorithm collects, from the tips to the
#' root, the set of states each node could have.  The second pass, from the
#' root to the tips, narrows each node's set to the states it shares with its
#' parent, if any.  Nodes left with more than one possible state are
#' ambiguous and get `NA`.  Missing tip values are treated as unknown (any
#' state), with a warning.  Branch lengths are ignored.
#'
#' @section Maximum likelihood:
#' The values are fitted with a continuous-time Markov model of state changes
#' along the branches, whose lengths are taken as time.  The `model` argument
#' is passed to [ape::ace()]: `"ER"` (one rate for all changes), `"SYM"`,
#' `"ARD"`, or a custom matrix.  In addition, `"ordered"` allows only changes
#' by one step between integer values, which suits chromosome numbers changed
#' by fusions and fissions; the states then include all the integers between
#' the smallest and the largest value, even if no species has them.
#'
#' Each internal node gets its most likely state if its marginal probability
#' is at least `threshold`, and `NA` otherwise.  Every tip needs a value.
#' The results are only as good as the tree: distance trees are not
#' time-calibrated phylogenies.
#'
#' @param tree A [`stbTree`] object.
#' @param values Either a named vector of values whose names are the tip labels
#'        of the tree, or a taxon table (a data frame with one row per species
#'        and the tip labels as row names).
#' @param column When `values` is a data frame, the name of the column
#'        containing the values.
#' @param method `"parsimony"` (default) or `"ML"` (maximum likelihood).
#' @param model For `method = "ML"`, the model of state changes (see below).
#' @param threshold For `method = "ML"`, the probability above which the most
#'        likely state of a node is reported.
#'
#' @returns The `tree` with a new column, named after `column` or `"state"`,
#' containing the tip values and the reconstructed internal node values, with
#' `NA` for ambiguous nodes.  With `method = "ML"`, a second column with the
#' `_prob` suffix contains the probability of the most likely state of each
#' node (1 for tips).
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Functions for trees
#' @seealso [visualizeTree()], which plots the reconstruction with `trait`.
#'
#' @examples
#' # Pretend chromosome numbers (the real ones are all 1).
#' chr <- c(Halobacterium_litoreum = 3, Halobacterium_noricense = 3,
#'          Halobacterium_salinarum = 3, Salarchaeum_japonicum  = 5,
#'          Haloferax_mediterranei  = 2, Haloferax_volcanii     = 2)
#' ancestralStates(Halo_Tree, chr)
#'
#' # With six species, maximum likelihood is uncertain everywhere.
#' ancestralStates(Halo_Tree, chr, method = "ML")
#'
#' @importFrom stats na.omit setNames
#' @export

ancestralStates <- function(tree, values, column = NULL, method = c("parsimony", "ML"),
                            model = "ER", threshold = 0.95) {
  method  <- match.arg(method)
  colname <- if (is.null(column)) "state" else column
  if (is.data.frame(values)) {
    if (is.null(column) || is.null(values[[column]]))
      stop("Column ", dQuote(column), " not found in the taxon table.")
    values <- setNames(values[[column]], rownames(values))
  }
  if (is.null(names(values)))
    stop("Values must be named after the tip labels of the tree.")

  tipLabels <- tree$label[tree$isTip]
  tipValues <- setNames(unname(values[tipLabels]), tipLabels)
  missing   <- tipLabels[is.na(tipValues)]

  if (method == "parsimony") {
    if (length(missing) > 0)
      warning("No value for: ", paste(missing, collapse = ", "), ".  Treating as unknown.")
    tree[[colname]] <- fitchStates(tree, tipValues)
  } else {
    if (length(missing) > 0)
      stop("Maximum likelihood needs a value for every tip.  Missing: ",
           paste(missing, collapse = ", "), ".")
    ml <- aceStates(tree, tipValues, model, threshold)
    tree[[colname]] <- ml$state
    tree[[paste0(colname, "_prob")]] <- ml$prob
  }
  tree
}

# Fitch parsimony.  Returns a vector of states in the order of the tree's rows.
fitchStates <- function(tree, tipValues) {
  tips  <- tree$node[tree$isTip]
  tipValues <- unname(tipValues[tree$label[tips]])
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
  result
}

# Maximum likelihood with ape::ace().  Returns a list of two vectors in the
# order of the tree's rows: the states (NA below threshold) and the
# probability of the most likely state.
#' @importFrom ape ace node.depth.edgelength
#' @importFrom tidytree as.phylo
#' @noRd
aceStates <- function(tree, tipValues, model, threshold) {
  phy <- as.phylo(tree)
  # Scale the tree to a height of 1.  This does not change the probabilities
  # (only the rates, which are not reported), but the optimiser of ape::ace()
  # can fail with very short or very long branches.
  phy$edge.length <- phy$edge.length / max(node.depth.edgelength(phy))
  x   <- tipValues[phy$tip.label]
  if (identical(model, "ordered")) {
    if (!is.numeric(x) || any(x != round(x)))
      stop("The ordered model needs integer values.")
    levels <- as.character(seq(min(x), max(x)))
    model  <- matrix(0L, length(levels), length(levels))
    model[abs(row(model) - col(model)) == 1] <- 1L
  } else {
    levels <- as.character(sort(unique(x)))
  }
  fit <- ace(factor(x, levels = levels), phy, type = "discrete", model = model) |>
    suppressWarnings() # ape warns when it can not compute standard errors.
  lik <- fit$lik.anc
  if (all(abs(lik - 1 / ncol(lik)) < 1e-3))
    warning("The model gave uniform probabilities to all states: it could ",
            "not be fitted to these data.")

  # Convert the state names back to the type of the input values.
  toType <- \(s) if (is.integer(tipValues)) as.integer(s)
                 else if (is.numeric(tipValues)) as.numeric(s) else s
  best <- toType(colnames(lik)[apply(lik, 1, which.max)])
  prob <- apply(lik, 1, max)
  best[prob < threshold] <- NA

  internal <- match(tree$node, as.integer(rownames(lik)))
  state <- best[internal]
  state[tree$isTip] <- unname(tipValues[tree$label[tree$isTip]])
  p <- prob[internal]
  p[tree$isTip] <- 1
  list(state = state, prob = unname(p))
}
