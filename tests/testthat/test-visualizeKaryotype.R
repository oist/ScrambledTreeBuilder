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

test_that("visualizeKaryotype() labels internal nodes with values", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), value = "Scrambling_index")
  i <- layer_index(p, "GeomLabel")
  d <- ggplot2::layer_data(p, i)
  expect_setequal(na.omit(d$label), round(na.omit(Halo_Tree$Scrambling_index), 2))
  # Same with a vector of values, and with another rounding.
  q <- visualizeKaryotype(Halo_Tree, halo_taxons(), value = Halo_Tree$Scrambling_index, valueround = 1)
  expect_setequal(na.omit(ggplot2::layer_data(q, i)$label), round(na.omit(Halo_Tree$Scrambling_index), 1))
})

test_that("visualizeKaryotype() node values have their own colour scale", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), value = "Scrambling_index")
  q <- visualizeKaryotype(Halo_Tree, halo_taxons())
  expect_equal(ggplot2::layer_data(p, 1)$colour, ggplot2::layer_data(q, 1)$colour)
  # ggnewscale renames the aesthetic of the karyotype scale.
  scales <- ggplot2::ggplot_build(p)$plot$scales$scales
  karyotype <- Filter(\(s) any(startsWith(s$aesthetics, "colour_ggnewscale")), scales)[[1]]
  expect_equal(karyotype$get_labels(), c("2", "3", "5", "ambiguous"))
  values <- Filter(\(s) identical(s$aesthetics, "colour"), scales)[[1]]
  expect_equal(values$name, "Scrambling_index")
  d <- ggplot2::layer_data(p, layer_index(p, "GeomLabel"))
  expect_gt(length(unique(d$colour[!is.na(d$label)])), 1)
})

test_that("visualizeKaryotype() can draw node values as points", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), value = "Scrambling_index", points = TRUE)
  expect_true(is.na(layer_index(p, "GeomLabel")))
  d <- ggplot2::layer_data(p, layer_index(p, "GeomPoint"))
  expect_equal(nrow(d), 5)  # Internal nodes only.
  expect_setequal(d$x, p$data$x[!p$data$isTip])
})

test_that("visualizeKaryotype() draws no node values by default", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  expect_true(is.na(layer_index(p, "GeomLabel")))
  expect_true(is.na(layer_index(p, "GeomPoint")))
})

test_that("visualizeKaryotype() axis shows pairwise distances from the tips", {
  # Drawn by default.
  p <- visualizeKaryotype(Halo_Tree, halo_taxons())
  x <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x
  xmax <- max(p$data$x[p$data$isTip])
  breaks <- x$get_breaks()
  labels <- as.numeric(x$get_labels())
  ok <- !is.na(breaks)
  expect_equal(labels[ok], 2 * (xmax - breaks[ok]))
  expect_true(0 %in% labels)
  # The Haloferax node sits at their (symmetrised) pairwise distance.
  node8 <- p$data$x[p$data$node == 8]
  expect_equal(2 * (xmax - node8), halo_sym()["Haloferax_mediterranei", "Haloferax_volcanii"])
  expect_s3_class(p$theme$axis.text.x, "ggplot2::element_text")
})

test_that("visualizeKaryotype() axis can show distances as percentages", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "percent")
  q <- visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "number")
  xp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x
  xq <- ggplot2::ggplot_build(q)$layout$panel_params[[1]]$x
  expect_equal(xp$get_breaks(), xq$get_breaks())
  ok <- !is.na(xq$get_breaks())
  expect_equal(xp$get_labels()[ok], paste0(100 * as.numeric(xq$get_labels()[ok]), "%"))
})

test_that("visualizeKaryotype() axis argument is partially matched", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "perc")
  expect_match(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels(), "%$",
               all = FALSE)
  expect_error(visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "yes"), "should be one of")
  expect_error(visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "n"), "should be one of")
  expect_error(visualizeKaryotype(Halo_Tree, halo_taxons(), axis = TRUE), "must be NULL or a character vector")
})

test_that("visualizeKaryotype() axis can be removed", {
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), axis = "none")
  expect_s3_class(p$theme$axis.text.x, "ggplot2::element_blank")
})

test_that("visualizeKaryotype() value labels have solid borders", {
  # The labels should not inherit the dashed linetype of ambiguous branches.
  p <- visualizeKaryotype(Halo_Tree, halo_taxons(), value = "Scrambling_index")
  d <- ggplot2::layer_data(p, layer_index(p, "GeomLabel"))
  expect_true(all(d$linetype %in% c("solid", 1)))
})
