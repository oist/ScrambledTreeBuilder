#' Add species traits to a tree
#'
#' Copy values that belong to each species, such as chromosome numbers, from
#' a taxon table to the tips of a tree, to plot them with the `trait` option
#' of [visualizeTree()].
#'
#' Species are matched by name, not by order.  Species of the table that are
#' not in the tree are ignored, so that one taxon table can serve several
#' trees.  Tips that are not in the table get `NA`, with a warning.  Internal
#' nodes always get `NA`: see [ancestralStates()] to reconstruct their values.
#'
#' @param tree A [`stbTree`] object.
#' @param traits Either a taxon table (a data frame with one row per species
#'        and the tip labels as row names), or a named vector of values whose
#'        names are the tip labels.
#' @param columns When `traits` is a data frame, the names of the columns to
#'        add.  By default, all of them.
#' @param colname When `traits` is a vector, the name of the new column.
#'
#' @returns The `tree` with new columns named after the traits, containing
#' their values on the tips and `NA` on the internal nodes.  Existing columns
#' with the same names are replaced.
#'
#' @author Claude Opus 5.5 (Anthropic)
#'
#' @family Functions for trees
#' @seealso [makeValueTibble()] to add values of pairs of species on the
#' internal nodes.
#'
#' @examples
#' Halo_Taxons
#' tree <- makeTraitTibble(Halo_Tree, Halo_Taxons)
#' visualizeTree(tree, value = NULL, trait = "Toy_trait")
#'
#' # A single trait from a named vector.
#' toy <- c(Haloferax_mediterranei = "x", Haloferax_volcanii = "y")
#' makeTraitTibble(subTree(Halo_Tree, Halo_FocalClades$Haloferax), toy, colname = "Toy_label")
#'
#' @export

makeTraitTibble <- function(tree, traits, columns = NULL, colname = "trait") {
  if (!is.data.frame(traits)) {
    if (is.null(names(traits)))
      stop("Values must be named after the tip labels of the tree.")
    traits <- setNames(data.frame(unname(traits), row.names = names(traits)), colname)
  }
  if (is.null(columns)) columns <- names(traits)
  absent <- setdiff(columns, names(traits))
  if (length(absent) > 0)
    stop("Column(s) not found in the taxon table: ", paste(dQuote(absent), collapse = ", "), ".")

  tips    <- tree$isTip
  rows    <- match(tree$label, rownames(traits))
  rows[!tips] <- NA
  missing <- tree$label[tips & is.na(rows)]
  if (length(missing) > 0)
    warning("No trait values for: ", paste(missing, collapse = ", "), ".")
  for (col in columns)
    tree[[col]] <- traits[[col]][rows]
  tree
}
