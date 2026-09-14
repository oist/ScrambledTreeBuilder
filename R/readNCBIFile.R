#' Gets the NCBI info for a given node.
#'
#' @param node Node
#' @param your_tibble Tibble data developed from tree.
#' @param your_ncbidf Dataframe containing the NCBI info from `ncbi_datasets.tsv`.
#' @param colname The desired column from the NCBI dataframe.
#'
#' @return The row in the dataframe corresponding to the specified node.
#'
#' @family Functions for trees
#'
#' @author Priscila Biller
#'
#' @export
#'
getNCBIinfo <- function(node, your_tibble, your_ncbidf, colname) {
  your_df[paste(gsub(" ", "_", your_ncbidf$Organism.Name), your_ncbidf$Assembly.Accession, sep = "_") == node,][[colname]]
}


#' Adds NCBI info to the tree
#'
#' Adds new node columns data in a tree, containing the NCBI information
#' saved in the file `ncbi_datasets.tsv`.
#'
#' @param your_tibble A [`ConvenientTblTree`] object.
#' @param input_file Path to the file `ncbi_datasets.tsv`.
#'
#' @return The `your_tibble` object with 3 additional columns:
#'   * `NCBI_id`, containing the NCBI identifier;
#'   * `NCBI_nbchr`, containing the number of chromosomes; and
#'   * `NCBI_seqlen`, containing the sequence length.
#'
#' @family Data load functions
#'
#' @author Priscila Biller
#'
#' @export
#'
readNCBIfile <- function(your_tibble, input_file) {
  # Read file with NCBI info.
  if (!file.exists(input_file)) {
    stop(paste("File does not exist: ", input_file))
  }
  df_ncbi <- read.table(file=input_file, sep='\t', header=TRUE)
  # Columns to be added for every node in the tree.
  desired_columns <- list(list(ncbi_name="Organism.Taxonomic.ID",                       tibble_name="NCBI_id"),
                          list(ncbi_name="Assembly.Stats.Total.Number.of.Chromosomes",  tibble_name="NCBI_nbchr"),
                          list(ncbi_name="Assembly.Stats.Total.Sequence.Length",        tibble_name="NCBI_seqlen"))
  # Add NCBI info for every node in the tree with a label (e.g. "Callithrix_jacchus_GCA_049354715.1")
  node_labels <- your_tibble$label[!is.na(your_tibble$label)] |> purrr::set_names()
  for(newcol in desired_columns){
    node_info <- node_labels |> sapply(getNCBIinfo, your_tibble, df_ncbi, newcol$ncbi_name)
    your_tibble[, newcol$tibble_name] <- NA
    your_tibble[match(names(node_info), your_tibble$label), newcol$tibble_name] <- unname(node_info)
  }
  your_tibble
}
