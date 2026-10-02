#' Pairwise comparison matrix
#'
#' Make a matrix from a data frame containing species data.
#'
#' @param pairwise_data Data frame containing species information, produced with [`formatStats()`]
#' @param column The variable from the data frame to extract species data from.
#' @param defaultDiagonal Maximum possible value of variable.
#' @param defaultValue Minimum possible value of variable.
#' @param impute Method to impute missing value. Defaults to `no`.
#' @param ... Other arguments passed to the imputing function.
#'
#' @note In case of missing values, the following options are available.
#'  - Leave them missing (`NA`). This is the default.
#'  - Give them a fixed value defined by the `devaultValue` option.
#'  - Impute them from the replicate (when _target_ and _query_ roles are
#'    exchanged) if available (`impute=average`), else give them the default
#'    value.
#'  - Impute them from the replicate, and then with random forests
#'    (`impute=missForest`).  Hide the diagonal values (which are extreme)
#'    when (`impute=missForest2`).
#'
#' @return A square `matrix` where rows and columns are sample identifiers and
#' values are statistics for the pairs of samples.  In addition it has an
#' attribute "builtWith" that records the value of `column`.
#'
#' @author Noa Brenner
#' @author Charles Plessy
#' @author Anika Mittal
#'
#' @importFrom matrixStats colAlls rowAlls
#'
#' @export
#'
#' @examples
#' makeMatrix(Halo_DF, "percent_difference_global", 0)
#'
#' # Missing values get NA by default unless specified in the 4th argument,
#' # or imputed from the reverse pair.
#' df <- Halo_DF
#' df["Salarchaeum_japonicum___Haloferax_volcanii", "percent_difference_global"] <- NA
#' makeMatrix(df, "percent_difference_global", 0, impute = "average")

makeMatrix <- function(pairwise_data, column="", defaultDiagonal = 100, defaultValue = NA, impute = c("no", "average", "missForest", "missForest2"), ...) {
  if(is.null(pairwise_data[[column]])) {
    warning("Column ", dQuote(column), " not found.")
    return(NULL)
  }
  impute <- match.arg(impute)
  m <- .makeMatrix(pairwise_data, column)
  problRow <- names(which(rowAlls(m == 0, na.rm = TRUE)))
  problCol <- names(which(colAlls(m == 0, na.rm = TRUE)))
  if (length(c(problRow, problCol)) !=0 )
    stop ("Problematic species: ", paste(unique(c(problRow, problCol)), collapse=", "))
  diag(m) <- defaultDiagonal
  m <- .imputeMatrix(m, defaultValue, impute, ...)
  attr(m, "builtWith") <- column
  m
}

.makeMatrix <- function(pairwise_data, column="") {
  all_species <- sort(unique(c(pairwise_data$species1, pairwise_data$species2)))
  m <- matrix(NA, nrow = length(all_species), ncol = length(all_species))
  colnames(m) <- rownames(m) <- all_species
  m[cbind(pairwise_data$species1, pairwise_data$species2)] <- pairwise_data[[column]]
  m
}

.imputeMatrix <- function(m, defaultValue, impute, ...) {

  fillSymmetricNA <- function(mat) {
    na_pos <- which(is.na(mat) & !is.na(t(mat)), arr.ind = TRUE)
    mat[na_pos] <- t(mat)[na_pos]
    mat
  }

  imputeDefault <- function(m, defaultValue, ...) {
    m[is.na(m)] <- defaultValue
    m
  }

  imputeAverage <- function(m, defaultValue, ...) {
    m <- fillSymmetricNA(m)
    imputeDefault(m, defaultValue, ...)
  }

  imputeForest  <- function(m, devaultValue, ...) {
    m <- fillSymmetricNA(m)
    if (!requireNamespace("missForest", quietly = TRUE))
      stop("Package 'missForest' is required for imputation. Please install it.")
    missForest::missForest(m, ...)$ximp |> suppressWarnings()
  }

  imputeForest2  <- function(m, devaultValue, ...) {
    diag.m.old <- diag(m)
    diag(m) <- NA
    m <- imputeForest(m, devaultValue, ...)
    diag(m) <- diag.m.old
    m
  }

  imputeFunction <- switch(impute,
                           no          = imputeDefault,
                           average     = imputeAverage,
                           missForest  = imputeForest,
                           missForest2 = imputeForest2)

  imputeFunction(m, defaultValue)
}
