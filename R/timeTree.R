#' Timeline - TimeTree
#' 
#' Fetches a list of divergence times of nodes in TimeTree from the specified taxon to the common ancestor of all cellular organisms. 
#' 
#' @param taxId An integer: The NCBI id of a taxon
#' @param divergenceTimeType Optional parameter specifying which type of divergence time to read from TimeTree. 
#'                    Two types are available; the default is `adjusted_age`:
#'                    1. `summary_age`: the median divergence time across all studies from TimeTree;
#'                    2. `adjusted_age`: a time adjusted when a node is older than its parent in the global timetree reconstructed from individual timetrees.
#' @param timetreeBaseTimeline Optional parameter specifying the base URL for retrieving timelines from the TimeTree API.
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
getTimeline <- function(taxId, divergenceTimeType="adjusted_age", timetreeBaseTimeline="http://timetree.temple.edu/api/timeline/") {
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
  timetreeRequest  <- paste0(timetreeBaseTimeline, taxId)
  timetreeResponse <- GET(timetreeRequest)
  lines  <- content(timetreeResponse, as="text", encoding="UTF-8")
  header <- scan(text=lines, what="", sep=",", nlines=1, quiet=TRUE) # First line: header; Other lines: data.
  # Check if the request returned something valid.
  if(length(header) < 8){return(NA)}
  # Parse request.
  df    <- read.csv(text=lines)
  # In TimeTree there are two possibilities for the divergence time (here we go with "adjusted age" as default):
  # - "Summary age" is the median of all studies.
  # - "Adjusted age" is the adjusted time, computed in cases where a node is older than its parent in the global timetree reconstructed from individual timetrees.
  df    <- df[ c("ncbi.id", "scientific_name", "rank", divergenceTimeType, "ci_low",         "ci_high")]
  names(df) <- c("NCBIid",  "name",            "rank", "divtime",           "divtimeCI_low", "divtimeCI_high")
  df
}

#' Divergence time - TimeTree
#' 
#' Retrieves the divergence time of two species via the TimeTree API.
#' 
#' @param taxId1 An integer: The NCBI id of one taxon
#' @param taxId2 An integer: The NCBI id of the other taxon
#' @param timetreeBasePairwise Optional parameter specifying the base URL for retrieving divergence times between pairs of species from the TimeTree API.
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
getDivergenceTime <- function(taxId1, taxId2, timetreeBasePairwise="http://timetree.temple.edu/api/pairwise/") {
  # Example of request to the TimeTree's API: 
  # - Request: http://timetree.temple.edu/api/pairwise/9555/9601
  # - Answer from TimeTree: 
  # taxon_a_id,taxon_b_id,scientific_name_a,scientific_name_b,all_total,precomputed_age,precomputed_ci_low,precomputed_ci_high,adjusted_age
  # 9555,9601,Papio anubis,Pongo abelii,83,28.82,26.8,30.6,0
  timetreeRequest  <- paste0(timetreeBasePairwise, taxId1, "/", taxId2)
  timetreeResponse <- GET(timetreeRequest)
  lines <- content(timetreeResponse, as="text", encoding="UTF-8")
  # Break the whole content into lines, and gets the second line (data).
  lineData <- strsplit(lines, "\n", fixed = TRUE)[[1]][2] # First line: header; Second line: data.
  # Break the line into data.
  data <- strsplit(lineData, ",", fixed = TRUE)[[1]]
  # Check if the request returned something valid.
  if(length(data) < 8){return(NA)}
  # Get data from request.
  list(divtime=as.numeric(data[6]),divtimeCI_low=as.numeric(data[7]),divtimeCI_high=as.numeric(data[8]))
}
