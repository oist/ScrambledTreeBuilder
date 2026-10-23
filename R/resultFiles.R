#' List YAML result files
#'
#' Given a directory containing the YAML files output by the pipeline
#' <https://github.com/oist/plessy_nf_GenomicBreaks/>, prepare a vector listing
#' these files, possibly removing the ones matching keywords give in a separate
#' input file.
#'
#' @note
#' The `remove` file is delimiter-separated, with the first column being a
#' pattern and the second an optional comment. It can be empty but must exist.
#'
#' @param dir Directory containing the YAML files.
#' @param remove File containing patterns for element removal.
#'
#' @returns A character vector where elements are file paths and names are
#' identifiers to be used by the `formatStats()` function.
#'
#' @author Charles Plessy
#' @author Priscila Biller
#'
#' @family Data load functions
#'
#' @examples
#' system.file("extdata/yaml", package = "ScrambledTreeBuilder") |> resultFiles()
#'
#' @importFrom utils read.delim
#' @export

resultFiles <- function (dir, remove = "removedAssemblies.txt") {
    files <- list.files(dir, pattern = ".yaml$|.yaml.bz2$|.yaml.gz$|.yaml.gz$", full.names = TRUE)
    names <- basename(files) |>
      sub(pattern = "\\.(gz|bz2|xz)$", replacement = "") |>
      sub(pattern = "\\.yaml",        replacement = "")
    names(files) <- names
    if(file.exists(remove)) {
      toRemove <- read.delim(remove, header = FALSE, col.names = c("pattern", "reason"))
      files    <- files[!grepl(paste(toRemove$pattern, collapse = "|"), files)]
    }
    files
}
