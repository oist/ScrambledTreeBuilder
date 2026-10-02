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

branch_colours <- function(p) setNames(ggplot2::layer_data(p, 1)$colour, p$data$node)

test_that("visualizeKaryotype() draws the most common state in grey", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  col <- branch_colours(p)
  expect_equal(col[["1"]], "grey30")                      # 3, on three tips
  expect_true(all(col[c("4", "6")] %in% grDevices::palette.colors(palette = "Okabe-Ito")))
})

test_that("visualizeKaryotype() draws ambiguous branches dashed and light grey", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  d <- ggplot2::layer_data(p, 1)
  ambiguous <- p$data$node %in% c(7, 9)
  expect_true(all(d$colour[ambiguous] == "grey70"))
  expect_true(all(d$linetype[ambiguous]  != "solid"))
  expect_true(all(d$linetype[!ambiguous] == "solid"))
})

test_that("visualizeKaryotype() background state and colors can be chosen", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), background = 2)
  col <- branch_colours(p)
  expect_equal(col[["4"]], "grey30")
  expect_false(col[["1"]] == "grey30")
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), colors = c(`2` = "red", `3` = "black", `5` = "blue"))
  col <- branch_colours(p)
  expect_equal(unname(col[c("1", "4", "6")]), c("black", "red", "blue"))
})

test_that("visualizeKaryotype() has enough colors for many states", {
  set.seed(1)
  m <- matrix(runif(100^2), 100, dimnames = list(paste0("sp", 1:100), paste0("sp", 1:100)))
  tree <- makeTidyTree((m + t(m)) / 2)
  taxons <- data.frame(row.names = paste0("sp", 1:100), ChromNumber = rep(1:20, each = 5))
  p <- visualizeKaryotype(tree, taxons)
  col <- branch_colours(p)[as.character(tree$node[tree$isTip])]
  expect_length(unique(col), 20)
})

test_that("visualizeKaryotype() passes reconstruction options to ancestralStates()", {
  tc <- two_clades()
  taxons <- data.frame(row.names = names(tc$values), ChromNumber = tc$values)
  p <- visualizeKaryotype(tc$tree, taxons, method = "ML", model = "ordered", threshold = 0.4)
  expect_true("ChromNumber_prob" %in% names(p$data))
  root <- p$data$node[p$data$parent == p$data$node]
  expect_equal(p$data$ChromNumber[p$data$node == root], 21)
})
