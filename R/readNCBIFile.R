#' Gets the NCBI info for a given node label.
#'
#' @param node_label The label from the desired node (e.g. "Callithrix_jacchus_GCA_049354715.1").
#' @param tree A [`ConvenientTblTree`] object.
#' @param ncbi_df Dataframe containing the NCBI info from `ncbi_datasets.tsv`.
#'
#' @return The row in the dataframe corresponding to the specified node label.
#'
#' @family Functions for trees
#'
#' @author Priscila Biller
#'
#' @export
#'
getNCBIinfo <- function(node_label, tree, ncbi_df) {
  ncbi_df[paste(gsub(" ", "_", ncbi_df$Organism.Name), ncbi_df$Assembly.Accession, sep = "_") == node_label,]
}


#' Adds NCBI info to the tree
#'
#' Adds new node columns data in a tree, containing the NCBI information
#' saved in the file `ncbi_datasets.tsv`.
#'
#' @param tree A [`ConvenientTblTree`] object.
#' @param ncbi_file Path to the file `ncbi_datasets.tsv`.
#'
#' @return The `tree` object with 3 additional columns:
#'   * `NCBI_id`, containing the NCBI identifier;
#'   * `NCBI_nbchr`, containing the number of chromosomes; and
#'   * `NCBI_seqlen`, containing the sequence length.
#'
#' @family Data load functions
#'
#' @author Priscila Biller
#'
#' @importFrom utils read.table
#' @export
#'
readNCBIfile <- function(tree, ncbi_file) {
  # Read file with NCBI info.
  if (!file.exists(ncbi_file)) {
    stop(paste("File does not exist: ", ncbi_file))
  }
  ncbi_df <- read.table(file=ncbi_file, sep='\t', header=TRUE)
  # Columns to be added for every node in the tree.
  desired_columns <- list(list(ncbi_name="Organism.Taxonomic.ID",                       tree_name="NCBI_id"),
                          list(ncbi_name="Assembly.Stats.Total.Number.of.Chromosomes",  tree_name="NCBI_nbchr"),
                          list(ncbi_name="Assembly.Stats.Total.Sequence.Length",        tree_name="NCBI_seqlen"))
  # Add NCBI info for every node in the tree with a label (e.g. "Callithrix_jacchus_GCA_049354715.1")
  nodes <- tree$label[!is.na(tree$label)] |> purrr::set_names() |> purrr::map(getNCBIinfo, tree=tree, ncbi_df=ncbi_df)
  print(nodes)
  for(newcol in desired_columns){
    tree[, newcol$tree_name] <- NA
    tree[match(names(nodes), tree$label), newcol$tree_name] <- (nodes |> purrr::map_vec(\(row) row[[newcol$ncbi_name]]))
  }
  tree
}
