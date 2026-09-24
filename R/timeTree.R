#' Timeline - TimeTree
#' 
#' Fetches a list of divergence times of nodes in TimeTree from the specified taxon to the common ancestor of all cellular organisms. 
#' 
#' @param taxid An integer: The NCBI id of a taxon
#'
#' @return A dataframe in which each row corresponds to an ancestor. 
#'         Each row contains the divergence time (`divtime`), its confidence interval 
#'         (`divtimeCI_low` and `divtimeCI_high`), and additional information, 
#'         including the scientific name (`name`) and taxonomic rank (`rank`).
#'         If the divergence time is not found in the TimeTree database, it returns `NA`.
#'
#' @author Priscila Biller
#'
#' @examples
#' \dontrun{
#' getTimeline(9555)}
#'
#' @importFrom httr GET content
#' @importFrom utils read.csv
#'
#' @export
getTimeline <- function(taxid) {
  timetreeBaseTimeline <- "http://timetree.temple.edu/api/timeline/"
  # Example of request to the TimeTree's API: 
  # - Request: http://timetree.temple.edu/api/timeline/9555
  # - Answer from TimeTree: 
  # ncbi-id,scientific_name,common_name,level,rank,summary_age,adjusted_age,ci_high,ci_low
  # -1,unnamed,,46,unknown,1.1,1.1,1.48,0.72
  # -1,unnamed,,45,unknown,1.27,1.36,1.97,0.56
  # -1,unnamed,,44,unknown,1.27,1.58,1.97,0.56
  # 9554,Papio,,43,genus,1.9,1.9,2.31,1.48
  # -1,unnamed,,42,unknown,4.68,4.68,5.08,4.27
  # ...
  timetreeRequest  <- paste(timetreeBaseTimeline, taxid, sep="")
  timetreeResponse <- GET(timetreeRequest)
  lines  <- content(timetreeResponse, as="text", encoding="UTF-8")
  header <- lines |> strsplit(split="\n", fixed=TRUE) |> (\(x) x[[1]][1])() |> strsplit(split = ",", fixed = TRUE) |> (\(x) x[[1]])() # First line: header; Other lines: data.
  if(length(header) >= 8) {
    df    <- read.csv(text=lines)
    # In TimeTree there are two possibilities for the divergence time (here we go with "adjusted age"):
    # - "Summary age" is the median of all studies.
    # - "Adjusted age" is the adjusted time, computed in cases where a node is older than its parent in the global timetree reconstructed from individual timetrees.
    df    <- df[ c("ncbi.id", "scientific_name", "rank", "adjusted_age", "ci_low",       "ci_high")]
    names(df) <- c("NCBIid",  "name",            "rank", "divtime",     "divtimeCI_low", "divtimeCI_high")
    df
  } else {
    NA
  }
}

#' Divergence time - TimeTree
#' 
#' Retrieves the divergence time of two species via the TimeTree API.
#' 
#' @param taxid1 An integer: The NCBI id of one taxon
#' @param taxid2 An integer: The NCBI id of the other taxon
#'
#' @return A list in which the first item is the divergence time (`divtime`),
#'         and the last two items are the divergence time confidence interval 
#'         (`divtimeCI_low` and `divtimeCI_high`).
#'         If the divergence time is not found in the TimeTree database, it returns `NA`.
#'
#' @author Priscila Biller
#'
#' @examples
#' \dontrun{
#' getDivergenceTime(9555,9601)}
#'
#' @importFrom httr GET content
#'
#' @export
getDivergenceTime <- function(taxid1, taxid2) {
  timetreeBasePairwise <- "http://timetree.temple.edu/api/pairwise/"
  # Example of request to the TimeTree's API: 
  # - Request: http://timetree.temple.edu/api/pairwise/9555/9601
  # - Answer from TimeTree: 
  # taxon_a_id,taxon_b_id,scientific_name_a,scientific_name_b,all_total,precomputed_age,precomputed_ci_low,precomputed_ci_high,adjusted_age
  # 9555,9601,Papio anubis,Pongo abelii,83,28.82,26.8,30.6,0
  timetreeRequest  <- paste(timetreeBasePairwise, taxid1, "/", taxid2, sep="")
  timetreeResponse <- GET(timetreeRequest)
  lines <- content(timetreeResponse, as="text", encoding="UTF-8")
  # Break the whole content into lines, and gets the second line (data).
  lineData <- strsplit(lines, "\n", fixed = TRUE)[[1]][2] # First line: header; Second line: data.
  # Break the line into data.
  data <- strsplit(lineData, ",", fixed = TRUE)[[1]]
  # Check if the request returned something valid.
  if(length(data) >= 8) {
    list(divtime=as.numeric(data[6]),divtimeCI_low=as.numeric(data[7]),divtimeCI_high=as.numeric(data[8]))
  } else {
    NA
  }
}
