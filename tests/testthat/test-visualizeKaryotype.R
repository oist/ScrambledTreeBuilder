# Index of the first layer whose geom is exactly of a class (ggtree's tip
# labels also inherit from GeomText).
layer_index <- function(p, class) which(sapply(p$layers, \(l) class(l$geom)[1] == class))[1]

test_that("visualizeKaryotype() colors branches by reconstructed state", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  expect_s3_class(p, "ggtree")
  expect_equal(p$data$ChromNumber, ancestralStates(Halo_Tree, halo_taxons(), "ChromNumber")$ChromNumber)
  # Horizontal branch segments, one per node.
  d <- ggplot2::layer_data(p, 1)
  col <- setNames(d$colour, p$data$node)
  expect_equal(unname(col[c("1", "2", "3", "10", "11")]), rep(col[["1"]], 5))  # all 3
  expect_equal(col[["4"]], col[["8"]])                                         # all 2
  expect_false(col[["4"]] == col[["1"]])
  expect_false(col[["6"]] %in% col[c("1", "4")])                               # 5
  expect_equal(unname(col[c("7", "9")]), rep(col[["7"]], 2))                   # ambiguous
  expect_false(col[["7"]] %in% col[c("1", "4", "6")])
})

test_that("visualizeKaryotype() prints the numbers next to the tips", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  d <- ggplot2::layer_data(p, layer_index(p, "GeomText"))
  expect_equal(nrow(d), 6)
  expect_equal(as.numeric(d$label[order(d$y)]),
               p$data$ChromNumber[p$data$isTip][order(p$data$y[p$data$isTip])])
  expect_true(all(d$x > max(p$data$x)))
  expect_length(unique(d$x), 1)
})

test_that("visualizeKaryotype() uses the column argument", {
  taxons <- halo_taxons()
  names(taxons)[2] <- "chr"
  p <- visualizeKaryotype(Halo_Tree, taxons, column = "chr")
  expect_true("chr" %in% names(p$data))
  expect_equal(p$labels$colour, "chr")
})

test_that("visualizeKaryotype() warns about species missing from the taxon table", {
  taxons <- halo_taxons()[-1, ]
  expect_warning(p <- visualizeKaryotype(Halo_Tree, taxons), "Halobacterium_litoreum")
  expect_s3_class(p, "ggtree")
})

test_that("visualizeKaryotype() plots can be highlighted with focal clades", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  expect_length((p + Halo_FocalClades)$layers, length(p$layers) + 2)
})

test_that("visualizeKaryotype() does not use deprecated ggplot2 arguments", {
  withr::local_options(lifecycle_verbosity = "error")
  expect_no_error(visualizeKaryotype(Halo_Tree, halo_taxons()))
})

test_that("visualizeKaryotype() labels ambiguous branches in the legend", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(scale$get_labels(), c("2", "3", "5", "ambiguous"))
})
