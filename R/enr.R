#' Compute ENR points (Excess Normalized Residuals over a LOESS trend)
#'
#' Trait-agnostic residual standardization for any y vs x relationship.
#' Fit a LOESS mean trend ŷ(x), compute residuals r = y - ŷ, and standardize
#' them to unit-free, z-like scores ENRz. Positive ENRz = above-trend; negative = below-trend.
#'
#' @param tb A data.frame/tibble with at least columns `x` and `y`.
#' @param span LOESS span for the mean trend (default 0.7).
#' @param family LOESS family; `"symmetric"` is robust to outliers (default).
#' @param degree Local polynomial degree for LOESS (1 = local lines).
#' @param use_weights If `TRUE`, use column `n` as LOESS weights when present.
#' @param std_mode One of `"robust_global"`, `"global_sd"`, or `"local_loess"`.
#'   - `robust_global`: global MAD*1.4826 (robust σ)
#'   - `global_sd`: global standard deviation
#'   - `local_loess`: LOESS on `abs(resid)` to model local σ
#' @param span_scale LOESS span used only for the `local_loess` sigma fit.
#'
#' @return `tb` with added columns:
#'   - `y_hat`: LOESS fitted value at each `x`
#'   - `resid`: raw residual `y - y_hat`
#'   - `ENRz`: standardized residual (z-like score)
#'
#' @examples
#' # MRCA summaries along a linear trend; clade A is shifted above it.
#' tb <- data.frame(x = 1:40, y = 0.02 * (1:40) + rnorm(40, sd = 0.02), n = 2,
#'                  clade = rep(c("A", "B", "C", "D"), each = 10))
#' tb$y[3:6] <- tb$y[3:6] + 0.15
#' computeENR(tb) |> head()
#'
#' @export
#' @importFrom stats loess predict mad sd
computeENR <- function(tb,
                       span = 0.7,
                       family = "symmetric",
                       degree = 1,
                       use_weights = TRUE,
                       std_mode = c("robust_global","global_sd","local_loess"),
                       span_scale = 0.7) {
  stopifnot(all(c("x","y") %in% names(tb)))
  std_mode <- match.arg(std_mode)
  
  # Optional LOESS weights
  w <- if (use_weights && "n" %in% names(tb)) tb$n else NULL
  
  # LOESS mean trend
  fit <- stats::loess(y ~ x, data = tb, span = span, family = family,
                      degree = degree, weights = w)
  
  # Fitted values and residuals
  tb$y_hat <- as.numeric(stats::predict(fit, newdata = data.frame(x = tb$x)))
  tb$resid <- tb$y - tb$y_hat
  
  # Standardization
  if (std_mode == "robust_global") {
    s <- stats::mad(tb$resid, center = 0, na.rm = TRUE) * 1.4826
    if (!is.finite(s) || s < .Machine$double.eps) s <- stats::sd(tb$resid, na.rm = TRUE)
    if (!is.finite(s) || s < .Machine$double.eps) s <- 1
    tb$ENRz <- tb$resid / s
    
  } else if (std_mode == "global_sd") {
    s <- stats::sd(tb$resid, na.rm = TRUE)
    if (!is.finite(s) || s < .Machine$double.eps) s <- 1
    tb$ENRz <- tb$resid / s
    
  } else { # local_loess
    abs_fit <- stats::loess(abs(resid) ~ x, data = tb, span = span_scale,
                            family = "symmetric")
    sigma_hat <- as.numeric(stats::predict(abs_fit, newdata = data.frame(x = tb$x))) * 1.253
    global_sd <- stats::sd(tb$resid, na.rm = TRUE); if (!is.finite(global_sd) || global_sd < .Machine$double.eps) global_sd <- 1
    bad <- !is.finite(sigma_hat) | sigma_hat < .Machine$double.eps
    sigma_hat[bad] <- global_sd
    tb$ENRz <- tb$resid / sigma_hat
  }
  
  tb
}

#' Aggregate ENR scores by clade
#'
#' Summarise point-level ENR statistics at the clade level using weighted means
#' (weights = `n` when available). Optionally require that "ABOVE-TREND" clades
#' also sit on the lower-divergence side (`x_mean` <= median `x_mean`).
#'
#' @param tb_points A data.frame/tibble containing columns `x,y,resid,ENRz,clade`
#'                  (and optionally `n`).
#' @param z_cut Threshold on standardized excess (`ENRz`) to tag categories (default 1).
#' @param require_low_x If `TRUE`, "ABOVE-TREND" must also have `x_mean` <= median `x_mean`.
#'
#' @return A tibble with one row per clade:
#'   `clade`, `ENR` (mean raw residual), `ENRz` (mean standardized),
#'   `x_mean`, `y_mean`, `n_nodes`, and `ENR_label`.
#'
#' @examples
#' tb <- data.frame(x = 1:40, y = 0.02 * (1:40) + rnorm(40, sd = 0.02), n = 2,
#'                  clade = rep(c("A", "B", "C", "D"), each = 10))
#' tb$y[3:6] <- tb$y[3:6] + 0.15
#' computeENR(tb) |> cladeENRtable(z_cut = 1)
#'
#' @export
#' @importFrom dplyr filter group_by summarise mutate case_when arrange desc n
#' @importFrom rlang .data
cladeENRtable <- function(tb_points, z_cut = 1, require_low_x = TRUE) {
  required <- c("x","y","resid","ENRz","clade")
  if (!all(required %in% names(tb_points))) {
    stop("tb_points must contain columns: x, y, resid, ENRz, clade")
  }
  
  clade_tbl <- tb_points |>
    dplyr::filter(!is.na(.data$clade), .data$clade != "Other") |>
    dplyr::group_by(.data$clade) |>
    dplyr::summarise(
      ENR     = stats::weighted.mean(.data$resid, .data$n, na.rm = TRUE),
      ENRz    = stats::weighted.mean(.data$ENRz, .data$n, na.rm = TRUE),
      x_mean  = stats::weighted.mean(.data$x, .data$n, na.rm = TRUE),
      y_mean  = stats::weighted.mean(.data$y, .data$n, na.rm = TRUE),
      n_nodes = sum(.data$n, na.rm = TRUE),
      .groups = "drop"
    )
  
  x_threshold <- stats::median(clade_tbl$x_mean, na.rm = TRUE)
  
  clade_tbl |>
    dplyr::mutate(
      ENR_label = dplyr::case_when(
        ENRz >=  z_cut & (!require_low_x | x_mean <= x_threshold) ~ "ABOVE-TREND",
        ENRz <= -z_cut                                            ~ "BELOW-TREND",
        TRUE                                                      ~ "TYPICAL"
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$ENRz))
}
