#' @include FocalClade.R
NULL

#' Subsample species evenly along evolutionary time
#'
#' Choose a subset of the species of a tree so that the most recent common
#' ancestors (MRCAs) that remain cover the range of a pairwise distance, used
#' as a proxy for evolutionary time, as evenly as possible, with replicates
#' from independent lineages at each time.  This reduces the number of pairs
#' of an all-versus-all comparison without losing time points.
#'
#' @section Method:
#' A subset of `n` species keeps exactly `n - 1` internal nodes of the tree,
#' and each species added to a subset adds exactly one node: the one where
#' it joins the species already kept.  The node values (for instance the
#' average pairwise distance between the species on each side of the node,
#' see [makeValueTibble()]) are cut in bins, and each bin gets a target
#' number of nodes: as equal as possible, but not more than the number of
#' nodes available in the bin.
#'
#' The subset starts with one species on each side of the root (or with the
#' species to `keep`), so that the deepest node is kept.  Species are then
#' added one by one, each time choosing the one whose new node falls in the
#' bin furthest below its target and, among these, in the focal clade least
#' represented in that bin, so that the replicates of a bin come from
#' independent lineages.  Finally, random swaps between kept and discarded
#' species are accepted when they reduce the shortfall to the targets, or the
#' number of nodes from the same clade in the same bin.
#'
#' The search is heuristic and uses random numbers: call [set.seed()] for
#' reproducible results.
#'
#' @section Choosing the number of species:
#' With `n`, the number of species is fixed and the `n - 1` nodes are spread
#' over the bins.  With `replicates`, the target of each bin is that number
#' of nodes (or all the nodes of the bin if it has fewer), and the number of
#' species is increased one by one until every target is met.  This lets the
#' shape of the tree decide how many species are needed: a bin whose nodes
#' are nested in each other costs more species than a bin whose nodes are in
#' separate lineages.
#'
#' @param tree A [`stbTree`] object.
#' @param value Name of a column of the tree containing the values of the
#'        internal nodes to bin, such as a pairwise distance projected with
#'        [makeValueTibble()], or a vector of these values.
#' @param n Number of species to keep.
#' @param replicates Number of nodes wanted in each bin, to choose the number
#'        of species automatically.  Give either `n` or `replicates`.
#' @param breaks Either the number of bins, of equal width from zero to the
#'        largest node value, or a vector of breaks covering the node values.
#' @param clades A [`FocalCladeList`] used to prefer nodes from different
#'        clades in the same bin.  Nodes outside the focal clades count as
#'        one clade.
#' @param keep Species that must be kept.
#' @param swaps Maximal number of swaps to evaluate after each addition of
#'        species.
#'
#' @returns A `SpeciesSubsample` object, with the properties `species` (the
#' kept species), `nodes` (their MRCAs), `report` (a data frame with the
#' number of nodes `available`, `target` and `kept` in each `bin`), and
#' `tree`, `value`, `valueName`, `breaks` and `clades` for
#' [plotSubsample()].
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Functions for trees
#' @seealso [plotSubsample()] to plot the result, [thinByMin()] to remove
#' near-identical genomes beforehand.
#'
#' @examples
#' set.seed(1)
#' (s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4, clades = Halo_FocalClades))
#' s@species
#'
#' # Let the tree decide how many species are needed for one node per bin.
#' subsampleSpecies(Halo_Tree, "Percent_difference", replicates = 1, breaks = 3)
#'
#' @export

subsampleSpecies <- function(tree, value, n = NULL, replicates = NULL, breaks = 10,
                             clades = NULL, keep = NULL, swaps = 20000) {
  if (is.null(n) == is.null(replicates))
    stop("Give either n or replicates.")
  valueName <- if (length(value) == 1 && is.character(value)) value else "value"
  if (length(value) == 1 && is.character(value)) {
    if (is.null(tree[[value]])) stop("Column ", dQuote(value), " not found in the tree.")
    value <- tree[[value]]
  }
  s <- subsampleSetup(tree, value, breaks, clades)
  nTips <- length(s$tips)
  if (!is.null(n) && (n < 2 || n > nTips))
    stop("n must be between 2 and ", nTips, ".")
  absent <- setdiff(keep, s$tips)
  if (length(absent) > 0)
    stop("Species to keep not found in the tree: ", paste(absent, collapse = ", "), ".")
  if (!is.null(n) && length(keep) > n)
    stop("More species to keep than n.")

  kept <- startingSpecies(s, keep)
  if (!is.null(n)) {
    target <- waterFill(s$available, n - 1)
    kept <- growSubsample(s, kept, max(n, length(kept)), target)
    kept <- swapSubsample(s, kept, keep, target, swaps)
  } else {
    target <- pmin(s$available, replicates)
    kept <- growSubsample(s, kept, max(sum(target) + 1, length(kept)), target)
    kept <- swapSubsample(s, kept, keep, target, swaps)
    while (sum(shortfall(s, kept, target)) > 0 && length(kept) < nTips) {
      kept <- growSubsample(s, kept, length(kept) + 1, target)
      kept <- swapSubsample(s, kept, keep, target, swaps)
    }
  }
  nodes <- keptNodes(s, kept)
  report <- data.frame(bin       = levels(s$bins),
                       available = as.vector(s$available),
                       target    = as.vector(target),
                       kept      = as.vector(table(s$bins[match(nodes, s$inner)])))
  SpeciesSubsample(species = kept, nodes = nodes, report = report, tree = tree,
                   value = value, valueName = valueName, breaks = s$breaks,
                   clades = clades)
}

# Result of subsampleSpecies().
SpeciesSubsample <- new_class("SpeciesSubsample", properties = list(
  species   = class_character,
  nodes     = class_numeric,
  report    = class_data.frame,
  tree      = class_any,
  value     = class_numeric,
  valueName = class_character,
  breaks    = class_numeric,
  clades    = class_any
))

method(print, SpeciesSubsample) <- function(x, ...) {
  cat(length(x@species), " of ", sum(x@tree$isTip), " species kept, with ",
      length(x@nodes), " MRCAs.  Nodes per bin of ", x@valueName, ":\n", sep = "")
  print(x@report, row.names = FALSE)
  missing <- sum(pmax(x@report$target - x@report$kept, 0))
  if (missing > 0) cat("Short of target by", missing, "node(s).\n")
  invisible(x)
}

# Precompute the tips under each side of each internal node, the bins, and
# the clade of each node.
subsampleSetup <- function(tree, value, breaks, clades) {
  tips  <- tree$label[tree$isTip]
  inner <- tree$node[!tree$isTip]
  sides <- lapply(inner, \(n) childSpecies(tree, n))
  L <- sapply(sides, \(x) tips %in% x$left)
  R <- sapply(sides, \(x) tips %in% x$right)
  dimnames(L) <- dimnames(R) <- list(tips, inner)
  time <- value[match(inner, tree$node)]
  if (anyNA(time)) stop("The values must be defined on all internal nodes.")
  if (length(breaks) == 1) breaks <- seq(0, max(time), length.out = breaks + 1)
  bins <- cut(time, breaks, include.lowest = TRUE)
  if (anyNA(bins)) stop("The breaks do not cover all the node values.")
  # Ancestors from the smallest to the largest.
  M <- (L | R)[, order(colSums(L | R)), drop = FALSE]
  nodeClade <- rep("Other", length(inner))
  if (!is.null(clades)) {
    for (cl in clades[order(-sapply(clades, \(cl) length(cl@genomeIDs)))]) {
      inside <- colSums((L | R) & !(tips %in% cl@genomeIDs)) == 0
      nodeClade[inside] <- cl@displayName
    }
  }
  list(tips = tips, inner = inner, L = L, R = R, M = M, bins = bins,
       breaks = breaks, available = table(bins), nodeClade = nodeClade)
}

# Spread a number of nodes over the bins as evenly as possible, without
# exceeding the number of nodes available in each bin.  When fewer nodes are
# left than open bins, they go to bins spread evenly from the deepest one, to
# keep the span of the values.
waterFill <- function(available, total) {
  target <- available * 0
  left <- total
  while (left > 0) {
    room <- available - target
    open <- which(room > 0)
    if (length(open) == 0) break
    if (left < length(open))
      open <- rev(open)[unique(round(seq(1, length(open), length.out = left)))]
    add <- pmin(room[open], max(1, left %/% length(open)))
    target[open] <- target[open] + add
    left <- left - sum(add)
  }
  target
}

# One species on each side of the root, unless the species to keep already
# span it.
startingSpecies <- function(s, keep) {
  root <- which.max(colSums(s$L | s$R))
  kept <- keep
  for (side in list(s$L[, root], s$R[, root]))
    if (!any(s$tips[side] %in% kept))
      kept <- c(kept, sample(s$tips[side], 1))
  kept
}

keptNodes <- function(s, kept) {
  k <- s$tips %in% kept
  s$inner[colSums(s$L & k) > 0 & colSums(s$R & k) > 0]
}

shortfall <- function(s, kept, target) {
  count <- table(s$bins[match(keptNodes(s, kept), s$inner)])
  pmax(target - count, 0)
}

subsampleObjective <- function(s, kept, target) {
  i <- match(keptNodes(s, kept), s$inner)
  same <- table(s$bins[i], s$nodeClade[i])
  1000 * sum(shortfall(s, kept, target)^2) + sum(choose(same, 2))
}

# Add species one by one, each time the one whose new node best fills the
# bins below target.  Looking one step ahead, a species also scores for the
# nodes between it and its join node: they need only one more species to be
# kept, which matters for nodes in lineages that have no kept species yet.
growSubsample <- function(s, kept, n, target) {
  ancestorNode <- match(colnames(s$M), s$inner)
  while (length(kept) < n) {
    cand    <- setdiff(s$tips, kept)
    hasKept <- colSums(s$M[kept, , drop = FALSE]) > 0
    onPath  <- s$M[cand, , drop = FALSE]
    joins   <- ancestorNode[apply(onPath & rep(hasKept, each = length(cand)), 1, \(x) which(x)[1])]
    i       <- match(keptNodes(s, kept), s$inner)
    count   <- table(s$bins[i])
    binDeficit <- (target - count) / pmax(target, 1)
    deficit <- binDeficit[s$bins[joins]]
    below   <- onPath & rep(!hasKept, each = length(cand))
    nodeDeficit <- rep(binDeficit[s$bins[ancestorNode]], each = length(cand))
    ahead   <- apply(ifelse(below, nodeDeficit, 0), 1, max)
    sameClade <- sapply(joins, \(j) sum(s$bins[i] == s$bins[j] & s$nodeClade[i] == s$nodeClade[j]))
    score   <- 100 * deficit + 50 * pmax(ahead, 0) - sameClade + stats::runif(length(cand), 0, 0.1)
    kept    <- c(kept, cand[which.max(score)])
  }
  kept
}

# Swap kept and discarded species while it improves the objective, in passes
# over the kept species in random order, until a pass brings no improvement
# or the number of evaluations reaches `swaps`.
swapSubsample <- function(s, kept, keep, target, swaps) {
  best  <- subsampleObjective(s, kept, target)
  evals <- 0
  repeat {
    improved <- FALSE
    for (i in sample(which(!kept %in% keep))) {
      for (sp in sample(setdiff(s$tips, kept))) {
        evals <- evals + 1
        if (evals > swaps) return(kept)
        k2 <- kept
        k2[i] <- sp
        o <- subsampleObjective(s, k2, target)
        if (o < best) {
          kept <- k2
          best <- o
          improved <- TRUE
          break
        }
      }
    }
    if (!improved) return(kept)
  }
}

#' Plot a subsample of species
#'
#' Show which species and which most recent common ancestors (MRCAs) are
#' kept by [subsampleSpecies()], how the MRCAs are spread over the bins of
#' node values, and, optionally, which pairs of an all-versus-all comparison
#' are kept.
#'
#' @param subsample A `SpeciesSubsample` object made by [subsampleSpecies()].
#' @param pairwise_data A data frame of pairwise results (one pair of
#'        genomes per row, from [formatStats()]), with `species1` and
#'        `species2` columns.
#' @param x,y Columns of `pairwise_data` to plot in the pairs panel.
#' @param panel `"all"` (default) to combine the panels with the
#'        _patchwork_ package, or the name of one panel: `"tree"` (kept
#'        species in red and kept MRCAs in gold), `"bins"` (histograms of
#'        the node values of all and kept MRCAs) or `"pairs"` (kept pairs in
#'        color, discarded pairs in grey).
#'
#' @returns A `ggplot` object, or a `patchwork` object for `panel = "all"`.
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Plotting functions
#' @seealso [subsampleSpecies()]
#'
#' @examples
#' set.seed(1)
#' s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4, clades = Halo_FocalClades)
#' plotSubsample(s, panel = "tree")
#' if (requireNamespace("patchwork", quietly = TRUE))
#'   plotSubsample(s, Halo_DF, x = "percent_difference_global")
#'
#' @importFrom ggplot2 aes facet_wrap geom_histogram geom_point labs
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual theme_bw
#' @importFrom ggtree geom_tippoint
#' @importFrom rlang .data
#' @export

plotSubsample <- function(subsample, pairwise_data = NULL, x = "percent_difference_local",
                          y = "index_avg_strandDiscord", panel = c("all", "tree", "bins", "pairs")) {
  panel <- match.arg(panel)
  s <- subsample
  if (panel == "pairs" && is.null(pairwise_data))
    stop("The pairs panel needs pairwise_data.")
  if (panel == "all") {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Combining the panels needs the patchwork package; or choose one panel.")
    right <- if (is.null(pairwise_data)) plotSubsample(s, panel = "bins")
             else patchwork::wrap_plots(plotSubsample(s, panel = "bins"),
                                        plotSubsample(s, pairwise_data, x, y, panel = "pairs"),
                                        ncol = 1)
    return(patchwork::wrap_plots(plotSubsample(s, panel = "tree"), right, widths = c(1.1, 1)))
  }

  cladeColors <- c(Other = "grey40")
  cladeOf <- \(sp) rep("Other", length(sp))
  if (!is.null(s@clades)) {
    for (cl in s@clades) cladeColors[cl@displayName] <- cl@color
    cladeOf <- function(sp) {
      out <- rep("Other", length(sp))
      for (cl in s@clades[order(-sapply(s@clades, \(cl) length(cl@genomeIDs)))])
        out[sp %in% cl@genomeIDs] <- cl@displayName
      out
    }
  }

  if (panel == "tree") {
    tree <- s@tree
    tree$Kept <- factor(ifelse(tree$isTip, ifelse(tree$label %in% s@species, "kept", "discarded"), NA),
                        c("kept", "discarded"))
    keptNode <- \(d) d[d$node %in% s@nodes, ]
    p <- visualizeTree(tree, value = NULL) +
      geom_tippoint(aes(colour = .data$Kept), size = 1.6) +
      geom_point(data = keptNode, aes(.data$x, .data$y), inherit.aes = FALSE,
                 fill = "gold", colour = "black", shape = 21, size = 2.2) +
      scale_colour_manual(values = c(kept = "firebrick", discarded = "grey75"),
                          name = "Species", na.translate = FALSE) +
      labs(title = sprintf("%d of %d species and %d MRCAs kept",
                           length(s@species), sum(tree$isTip), length(s@nodes)))
    if (!is.null(s@clades)) p <- p + cladeBars(s@clades)
    return(p)
  }

  if (panel == "bins") {
    inner <- s@tree$node[!s@tree$isTip]
    nodeClade <- vapply(inner, \(n) {
      sp <- unlist(childSpecies(s@tree, n))
      cl <- unique(cladeOf(sp))
      if (length(cl) == 1) cl else "Other"
    }, character(1))
    all  <- data.frame(value = s@value[match(inner, s@tree$node)], clade = nodeClade)
    sets <- c("all MRCAs", "kept MRCAs")
    d <- rbind(cbind(all, set = sets[1]), cbind(all[inner %in% s@nodes, ], set = sets[2]))
    d$set <- factor(d$set, sets)
    return(ggplot2::ggplot(d, aes(.data$value, fill = .data$clade)) +
      geom_histogram(breaks = s@breaks, colour = "white") +
      facet_wrap(~ set, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = cladeColors, name = "Clade of MRCA") +
      labs(title = "MRCAs per bin", x = s@valueName, y = "MRCAs") +
      theme_bw())
  }

  needed <- c("species1", "species2", x, y)
  absent <- setdiff(needed, names(pairwise_data))
  if (length(absent) > 0)
    stop("Column(s) not found in pairwise_data: ", paste(absent, collapse = ", "), ".")
  d <- pairwise_data[, needed]
  d$kept  <- d$species1 %in% s@species & d$species2 %in% s@species
  c1 <- cladeOf(d$species1)
  d$clade <- ifelse(c1 == cladeOf(d$species2), c1, "Other")
  ggplot2::ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_point(data = d[!d$kept, ], colour = "grey85", size = 0.6) +
    geom_point(data = d[d$kept, ], aes(colour = .data$clade), size = 1.1, alpha = 0.8) +
    scale_colour_manual(values = cladeColors, name = "Pair within") +
    labs(title = sprintf("%d of %d pairs kept", sum(d$kept), nrow(d)), x = x, y = y) +
    theme_bw()
}
