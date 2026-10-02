#' Most recent common ancestors
#'
#' Summarise statistics for MRCAs
#'
#' @param pairwise_data A results data frame produced by [`formatStats()`], with MRCAs
#'        annotated with [`recordAncestor()`].
#' @param clades A [`FocalCladeList`] object.
#' @param x,y the name of the pairwise statistics to summarise.
#' @param center,dispersion Functions to compute a representation of the
#'        central tendency and the dispersion of the values for each MRCA.
#'
#' @author Charles Plessy, Takahiro Fujita
#'
#' @returns A tibble ([tibble::tbl_df-class]), the raw and aggregated
#' information provided as one pair of `x` and `y` coordinates per line.
#'
#' This function is run directly by [`MRCA_2D_plot`], therefore it is not
#' exported anymore.
#'
#' @importFrom dplyr group_by summarize sym ungroup filter transmute slice_max left_join mutate
#' @importFrom tidyselect all_of
#' @importFrom tibble tibble

MRCAs <- function(pairwise_data, clades = NULL, x = "percent_difference_local", y = "index_avg_strandDiscord", center = mean, dispersion = sd) {

  # TODO: Use stbPairwiseData S7 class for validation (see issue #XXX)
  stopifnot("MRCA" %in% names(pairwise_data))

  tb <- pairwise_data |>
    group_by(.data$MRCA) |>
    summarize(
      x    = center(     !!sym(x), na.rm = TRUE),
      y    = center(     !!sym(y), na.rm = TRUE),
      xerr = dispersion( !!sym(x), na.rm = TRUE),
      yerr = dispersion( !!sym(y), na.rm = TRUE),
      n = length(.data$MRCA)
    ) |>
    ungroup()

  tb$clade <- "Other"
  tb$color <- NA
  if (! is.null(clades)) {
    for (clade in clades) {
      tb$clade[tb$MRCA %in% clade@nodeList] <- sub(" .*", "", clade@displayName)
      tb$color[tb$MRCA %in% clade@nodeList] <- sub(" .*", "", clade@color)
    }
  }

  tb$hover_text <- paste0(
    "MRCA: ", tb$MRCA,
    "<br>Clade: ", tb$clade,
    "<br>x: ", round(tb$x, 2),
    "<br>y: ", round(tb$y, 2),
    "<br>n: ", tb$n
  )
  tb$type <- "MRCA"

  # if (!is.null(mrca_label)) {
  #   stopifnot(mrca_label == "representative_pair")
  #   stopifnot(y %in% names(pairwise_data))
  #   stopifnot("species1" %in% names(pairwise_data))
  #   stopifnot("species2" %in% names(pairwise_data))
  #
  #   rep_pair <- pairwise_data |>
  #     filter(!is.na(.data$MRCA),
  #            !is.na(.data[[y]]),
  #            !is.na(.data$species1),
  #            !is.na(.data$species2)) |>
  #     group_by(.data$MRCA) |>
  #     slice_max(order_by = .data[[y]], n = 1, with_ties = FALSE) |>
  #     ungroup() |>
  #     transmute(
  #       MRCA = .data$MRCA,
  #       mrca_label = paste0(.data$species1, " vs ", .data$species2)
  #     )
  #
  #   tb <- tb |> left_join(rep_pair, by = "MRCA")
  # }

  tb <- tb[tb$MRCA != 0,]

  pairs <- tibble(
    MRCA  = pairwise_data$MRCA,
    x     = pairwise_data[[x]],
    y     = pairwise_data[[y]],
    xerr  = NA,
    yerr  = NA,
    n     = NA,
    clade = pairwise_data$focalClade,
    color = pairwise_data$focalColor,
    hover_text = paste0(
      "Pair: ", pairwise_data$species1, " vs ", pairwise_data$species2,
      "<br>x: ", round(pairwise_data[[.env$x]], 2),
      "<br>y: ", round(pairwise_data[[.env$y]], 2),
      "<br>Clade: ", pairwise_data$focalClade),
    type  = "Pair"
  )
  pairs[is.na(pairs$clade), "clade"] <- "Other"

  rbind(tb, pairs)
}

#' 2D MRCA plot
#'
#' 2D plot of the statistics extracted for each MRCA.  Optionally plots error
#' bars and the pair data in addition to its aggregated statistics.
#'
#' @param pairwise_data A results data frame produced by [`formatStats()`], with MRCAs
#'        annotated with [`recordAncestor()`].
#' @param clades A [`FocalCladeList`] object.
#' @param x,y Column names for the statistics to plot.
#' @param errorbars Plot the error bars?
#' @param pairs Plot the pairs?
#' @param xlim,ylim Maximal X and Y value.
#'
#' @returns A [`ggplot2::ggplot2`] object in which some information was added
#' for the [`plotly::ggplotly`] function to generate hover information.
#'
#' @author Charles Plessy, Takahiro Fujita
#'
#' @family Plotting functions
#' @seealso [averageResults()]
#'
#' @importFrom ggplot2 aes ggplot geom_errorbar geom_errorbarh geom_point labs
#' @importFrom ggplot2 scale_color_manual scale_x_continuous scale_y_continuous theme_minimal
#' @importFrom ggplot2 scale_size_identity scale_alpha_identity
#' @importFrom plotly ggplotly
#' @importFrom rlang .data .env
#' @importFrom stats na.omit sd
#' @export
#'
#' @examples
#' MRCA_2D_plot(Halo_DF, Halo_FocalClades) + ggplot2::labs(x="nucl. diff", y="scrambling")
#' MRCA_2D_plot(averageResults(Halo_DF), Halo_FocalClades)
#' MRCA_2D_plot(Halo_DF, Halo_FocalClades) |> plotly::ggplotly(tooltip ="text")

MRCA_2D_plot <- function(pairwise_data, clades = NULL, x = "percent_difference_local", y = "index_avg_strandDiscord", errorbars = TRUE, pairs = TRUE, xlim = 40, ylim = 1) {
  plot_data       <- MRCAs(pairwise_data, clades, x = x, y = y)
  plot_data$size  <- ifelse(plot_data$type == "MRCA", 4, 2.0)
  plot_data$alpha <- ifelse(plot_data$type == "MRCA", 1, 0.7)

  if (! isTRUE(pairs)) plot_data <- plot_data[plot_data$type == "MRCA",]

  p <- ggplot(plot_data,
              aes(x    = .data$x,
                  y    = .data$y,
                  col  = .data$clade,
                  text = .data$hover_text)) +
         geom_point(aes(size=.data$size, alpha = .data$alpha)) +
         scale_size_identity() +
         scale_alpha_identity() +
         labs(subtitle = "Average statistics for most recent common ancestors") +
    theme_minimal() +
    scale_x_continuous(limits = c(0, xlim)) +
    scale_y_continuous(limits = c(0, ylim))

  ## Manual color scale (same logic you already had)
  colors.tb <- plot_data |>
    dplyr::group_by(.data$clade) |>
    dplyr::summarize(
      clade = unique(.data$clade),
      color = unique(.data$color),
      .groups = "drop"
    )
  colors <- colors.tb$color
  names(colors) <- colors.tb$clade
  colors <- stats::na.omit(colors)
  p <- p + scale_color_manual(values = colors)

  ## Error bars ONLY for MRCAs
  if (isTRUE(errorbars)) {
    p <- p +
      geom_errorbar(
        data = dplyr::filter(plot_data, .data$type == "MRCA"),
        aes(ymin = .data$y - .data$yerr,
            ymax = .data$y + .data$yerr),
        width = 0
      ) +
      geom_errorbar(
        data = dplyr::filter(plot_data, .data$type == "MRCA"),
        aes(xmin = .data$x - .data$xerr,
            xmax = .data$x + .data$xerr,
            y = .data$y),
        orientation = "y",
        width = 0
      )
  }

  p
}
