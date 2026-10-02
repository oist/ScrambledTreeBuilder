test_that("MRCAs() summarises pairs by most recent common ancestor", {
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF, Halo_FocalClades)
  expect_named(tb, c("MRCA", "x", "y", "xerr", "yerr", "n", "clade", "color", "hover_text", "type"))
  mrca <- tb[tb$type == "MRCA", ]
  expect_equal(mrca$MRCA, c(7, 8, 9, 10, 11))
  expect_equal(mrca$n,    c(16L, 2L, 6L, 4L, 2L))
  pairs8 <- Halo_DF[Halo_DF$MRCA == 8, ]
  expect_equal(mrca$x[mrca$MRCA == 8], mean(pairs8$percent_difference_local))
  expect_equal(mrca$y[mrca$MRCA == 8], mean(pairs8$index_avg_strandDiscord))
  expect_equal(mrca$xerr[mrca$MRCA == 8], sd(pairs8$percent_difference_local))
})

test_that("MRCAs() adds one row per pair", {
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF, Halo_FocalClades)
  pairs <- tb[tb$type == "Pair", ]
  expect_equal(nrow(pairs), nrow(Halo_DF))
  expect_equal(pairs$x, Halo_DF$percent_difference_local)
  expect_true(all(is.na(pairs$xerr)))
  expect_true(all(startsWith(pairs$hover_text, "Pair: ")))
})

test_that("MRCAs() colors MRCAs by focal clade", {
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF, Halo_FocalClades)
  mrca <- tb[tb$type == "MRCA", ]
  expect_equal(mrca$clade, c("Other", "Haloferax", "Other", "Halobacterium", "Halobacterium"))
  expect_equal(mrca$color, c(NA, "green3", NA, "blue", "blue"))
  expect_equal(mrca$hover_text[2],
               paste0("MRCA: 8<br>Clade: Haloferax<br>x: ", round(mrca$x[2], 2),
                      "<br>y: ", round(mrca$y[2], 2), "<br>n: 2"))
})

test_that("MRCAs() works without focal clades", {
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF)
  expect_true(all(tb$clade[tb$type == "MRCA"] == "Other"))
})

test_that("MRCAs() uses the requested statistics and summary functions", {
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF, x = "percent_difference_global",
                                     y = "percent_aligned", center = median)
  pairs7 <- Halo_DF[Halo_DF$MRCA == 7, ]
  expect_equal(tb$x[tb$type == "MRCA" & tb$MRCA == 7], median(pairs7$percent_difference_global))
  expect_equal(tb$y[tb$type == "MRCA" & tb$MRCA == 7], median(pairs7$percent_aligned))
})

test_that("MRCAs() drops the aggregate of pairs absent from the tree", {
  df <- recordAncestor(Halo_DF, subTree(Halo_Tree, 9))
  tb <- ScrambledTreeBuilder:::MRCAs(df)
  expect_false(0 %in% tb$MRCA[tb$type == "MRCA"])
})

test_that("MRCAs() requires an MRCA column", {
  expect_error(ScrambledTreeBuilder:::MRCAs(formatStats(halo_files())), "MRCA")
})

test_that("MRCA_2D_plot() plots MRCAs, pairs and error bars", {
  p <- MRCA_2D_plot(Halo_DF, Halo_FocalClades)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 35)
  expect_length(p$layers, 3)
  expect_equal(p$data$size[p$data$type == "MRCA"], rep(4, 5))
  expect_equal(p$data$alpha[p$data$type == "Pair"], rep(0.7, 30))
})

test_that("MRCA_2D_plot() options remove pairs and error bars", {
  expect_equal(nrow(MRCA_2D_plot(Halo_DF, Halo_FocalClades, pairs = FALSE)$data), 5)
  expect_length(MRCA_2D_plot(Halo_DF, Halo_FocalClades, errorbars = FALSE)$layers, 1)
})

test_that("MRCA_2D_plot() passes the statistics to plot", {
  p <- MRCA_2D_plot(Halo_DF, x = "percent_difference_global", y = "percent_aligned", xlim = 100, ylim = 100)
  pairs <- p$data[p$data$type == "Pair", ]
  expect_equal(pairs$x, Halo_DF$percent_difference_global)
  expect_equal(pairs$y, Halo_DF$percent_aligned)
})

test_that("MRCA_2D_plot() uses focal clade colors", {
  p <- MRCA_2D_plot(Halo_DF, Halo_FocalClades)
  b <- ggplot2::layer_data(p, 1)
  expect_true(all(b$colour[p$data$clade == "Haloferax"] == "green3"))
  expect_true(all(b$colour[p$data$clade == "Halobacterium"] == "blue"))
})

test_that("MRCA_2D_plot() output can be converted by plotly", {
  expect_s3_class(plotly::ggplotly(MRCA_2D_plot(Halo_DF, Halo_FocalClades), tooltip = "text"), "plotly")
})
