# Index of the first layer whose geom is exactly of a class (ggtree's tip
# labels also inherit from GeomText).
layer_index <- function(p, class) which(sapply(p$layers, \(l) class(l$geom)[1] == class))[1]

test_that("visualizeTree() labels internal nodes with their IDs by default", {
  p <- visualizeTree(Halo_Tree)
  expect_s3_class(p, "ggtree")
  b <- ggplot2::ggplot_build(p)
  labels <- b$data[[length(b$data)]]$label
  expect_setequal(labels, Halo_Tree$node)
  expect_null(ggplot2::get_guide_data(p, "colour"))
})

test_that("visualizeTree() hides the node ID legend but not the trait legend", {
  p <- visualizeTree(halo_chr_tree(), trait = "ChromNumber")
  expect_setequal(na.omit(ggplot2::layer_data(p, layer_index(p, "GeomLabel"))$label), Halo_Tree$node)
  expect_null(ggplot2::get_guide_data(p, "colour"))
  expect_equal(ggplot2::get_guide_data(p, "colour_ggnewscale_1")$.label, c("2", "3", "5", "ambiguous"))
})

test_that("visualizeTree() names the legend after the value column", {
  p <- visualizeTree(Halo_Tree, "Scrambling_index")
  expect_equal(p$labels$colour, "Scrambling_index")
})

test_that("visualizeTree() draws nothing on nodes with value = NULL", {
  p <- visualizeTree(Halo_Tree, value = NULL)
  expect_true(is.na(layer_index(p, "GeomLabel")))
})

test_that("visualizeTree() labels nodes with a column or a vector of values", {
  by_name  <- visualizeTree(Halo_Tree, "Scrambling_index")
  by_value <- visualizeTree(Halo_Tree, Halo_Tree$Scrambling_index)
  b <- ggplot2::ggplot_build(by_name)
  labels <- b$data[[length(b$data)]]$label
  expect_setequal(na.omit(labels), round(na.omit(Halo_Tree$Scrambling_index), 2))
  expect_equal(ggplot2::layer_data(by_name, length(by_name$layers))$label,
               ggplot2::layer_data(by_value, length(by_value$layers))$label)
})

test_that("visualizeTree() rounds and nudges labels", {
  p <- visualizeTree(Halo_Tree, "Percent_difference", valueround = 0, ynudge = 0.5)
  q <- visualizeTree(Halo_Tree, "Percent_difference", valueround = 0)
  n <- length(p$layers)
  labels <- ggplot2::layer_data(p, n)$label
  expect_true(all(na.omit(labels) == round(na.omit(labels))))
  expect_equal(ggplot2::layer_data(p, n)$y - ggplot2::layer_data(q, n)$y, rep(0.5, 11))
})

test_that("visualizeTree() sets the label border width", {
  p <- visualizeTree(Halo_Tree, outerlabelsize = 0.7)
  expect_equal(p$layers[[length(p$layers)]]$aes_params$linewidth, 0.7)
})

test_that("visualizeTree() does not use deprecated ggplot2 arguments", {
  withr::local_options(lifecycle_verbosity = "error")
  expect_no_error(visualizeTree(Halo_Tree))
})

test_that("visualizeTree() output is unchanged by the internal helper", {
  p <- visualizeTree(Halo_Tree, "Scrambling_index", ynudge = 0.2)
  l <- p$layers[[length(p$layers)]]
  expect_s3_class(l$geom, "GeomLabel")
  expect_equal(l$aes_params$linewidth, 0.25)
  expect_equal(l$aes_params$size, 3)
})

test_that("visualizeTree() colors branches by reconstructed state", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
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

test_that("visualizeTree() prints the numbers next to the tips", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  d <- ggplot2::layer_data(p, layer_index(p, "GeomText"))
  expect_equal(nrow(d), 6)
  expect_equal(as.numeric(d$label[order(d$y)]),
               p$data$ChromNumber[p$data$isTip][order(p$data$y[p$data$isTip])])
  expect_true(all(d$x > max(p$data$x)))
  expect_length(unique(d$x), 1)
})

test_that("visualizeTree() trait can be any column of the tree", {
  tree <- halo_chr_tree()
  names(tree)[names(tree) == "ChromNumber"] <- "chr"
  p <- visualizeTree(tree, value = NULL, trait = "chr")
  expect_true("chr" %in% names(p$data))
  expect_equal(p$labels$colour, "chr")
})

test_that("visualizeTree() trait ignores values on internal nodes", {
  tree <- halo_chr_tree()
  tree$ChromNumber[!tree$isTip] <- 99
  expect_equal(visualizeTree(tree, trait = "ChromNumber")$data$ChromNumber,
               visualizeTree(halo_chr_tree(), trait = "ChromNumber")$data$ChromNumber)
})

test_that("visualizeTree() warns about tips without trait value", {
  tree <- halo_chr_tree()
  tree$ChromNumber[1] <- NA
  expect_warning(p <- visualizeTree(tree, trait = "ChromNumber"), "Halobacterium_litoreum")
  expect_s3_class(p, "ggtree")
})

test_that("visualizeTree() trait must be a column of the tree", {
  expect_error(visualizeTree(Halo_Tree, trait = "ChromNumber"), "ChromNumber")
  expect_error(visualizeTree(halo_chr_tree(), trait = c("ChromNumber", "label")), "name of a column")
  expect_error(visualizeTree(halo_chr_tree(), trait = halo_taxons()), "name of a column")
})

test_that("visualizeTree() plots can be highlighted with focal clades", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  expect_length((p + Halo_FocalClades)$layers, length(p$layers) + 2)
})

test_that("visualizeTree() trait plots do not use deprecated ggplot2 arguments", {
  withr::local_options(lifecycle_verbosity = "error")
  expect_no_error(visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber"))
})

test_that("visualizeTree() labels ambiguous branches in the legend", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(scale$get_labels(), c("2", "3", "5", "ambiguous"))
})

branch_colours <- function(p) setNames(ggplot2::layer_data(p, 1)$colour, p$data$node)

test_that("visualizeTree() draws the most common state in grey", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  col <- branch_colours(p)
  expect_equal(col[["1"]], "grey30")                      # 3, on three tips
  expect_true(all(col[c("4", "6")] %in% grDevices::palette.colors(palette = "Okabe-Ito")))
})

test_that("visualizeTree() draws ambiguous branches dashed and light grey", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  d <- ggplot2::layer_data(p, 1)
  ambiguous <- p$data$node %in% c(7, 9)
  expect_true(all(d$colour[ambiguous] == "grey70"))
  expect_true(all(d$linetype[ambiguous]  != "solid"))
  expect_true(all(d$linetype[!ambiguous] == "solid"))
})

test_that("visualizeTree() background state and colors can be chosen", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber", background = 2)
  col <- branch_colours(p)
  expect_equal(col[["4"]], "grey30")
  expect_false(col[["1"]] == "grey30")
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber", colors = c(`2` = "red", `3` = "black", `5` = "blue"))
  col <- branch_colours(p)
  expect_equal(unname(col[c("1", "4", "6")]), c("black", "red", "blue"))
})

test_that("visualizeTree() has enough colors for many states", {
  set.seed(1)
  m <- matrix(runif(100^2), 100, dimnames = list(paste0("sp", 1:100), paste0("sp", 1:100)))
  tree <- makeTidyTree((m + t(m)) / 2)
  tree$ChromNumber <- setNames(rep(1:20, each = 5), paste0("sp", 1:100))[tree$label]
  p <- visualizeTree(tree, value = NULL, trait = "ChromNumber")
  col <- branch_colours(p)[as.character(tree$node[tree$isTip])]
  expect_length(unique(col), 20)
})

test_that("visualizeTree() passes reconstruction options to ancestralStates()", {
  tc <- two_clades()
  tree <- tc$tree
  tree$ChromNumber <- tc$values[tree$label]
  p <- visualizeTree(tree, trait = "ChromNumber", method = "ML", model = "ordered", threshold = 0.4)
  expect_true("ChromNumber_prob" %in% names(p$data))
  root <- p$data$node[p$data$parent == p$data$node]
  expect_equal(p$data$ChromNumber[p$data$node == root], 21)
})

test_that("visualizeTree() labels internal nodes with values", {
  p <- visualizeTree(halo_chr_tree(), "Scrambling_index", trait = "ChromNumber")
  i <- layer_index(p, "GeomLabel")
  d <- ggplot2::layer_data(p, i)
  expect_setequal(na.omit(d$label), round(na.omit(Halo_Tree$Scrambling_index), 2))
  # Same with a vector of values, and with another rounding.
  q <- visualizeTree(halo_chr_tree(), Halo_Tree$Scrambling_index, trait = "ChromNumber", valueround = 1)
  expect_setequal(na.omit(ggplot2::layer_data(q, i)$label), round(na.omit(Halo_Tree$Scrambling_index), 1))
})

test_that("visualizeTree() node values have their own colour scale", {
  p <- visualizeTree(halo_chr_tree(), "Scrambling_index", trait = "ChromNumber")
  q <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  expect_equal(ggplot2::layer_data(p, 1)$colour, ggplot2::layer_data(q, 1)$colour)
  # ggnewscale renames the aesthetic of the trait scale.
  scales <- ggplot2::ggplot_build(p)$plot$scales$scales
  trait <- Filter(\(s) any(startsWith(s$aesthetics, "colour_ggnewscale")), scales)[[1]]
  expect_equal(trait$get_labels(), c("2", "3", "5", "ambiguous"))
  expect_equal(ggplot2::get_labs(p)$colour, "Scrambling_index")
  expect_equal(ggplot2::get_labs(p)[[trait$aesthetics[1]]], "ChromNumber", ignore_attr = TRUE)
  d <- ggplot2::layer_data(p, layer_index(p, "GeomLabel"))
  expect_gt(length(unique(d$colour[!is.na(d$label)])), 1)
})

test_that("visualizeTree() can draw node values as points", {
  p <- visualizeTree(halo_chr_tree(), "Scrambling_index", trait = "ChromNumber", points = TRUE)
  expect_true(is.na(layer_index(p, "GeomLabel")))
  d <- ggplot2::layer_data(p, layer_index(p, "GeomPoint"))
  expect_equal(nrow(d), 5)  # Internal nodes only.
  expect_setequal(d$x, p$data$x[!p$data$isTip])
})

test_that("visualizeTree() draws no node values with value = NULL", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber")
  expect_true(is.na(layer_index(p, "GeomLabel")))
  expect_true(is.na(layer_index(p, "GeomPoint")))
})

test_that("visualizeTree() axis shows pairwise distances from the tips", {
  # Drawn by default.
  p <- visualizeTree(Halo_Tree)
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

test_that("visualizeTree() axis can show distances as percentages", {
  p <- visualizeTree(Halo_Tree, axis = "percent")
  q <- visualizeTree(Halo_Tree, axis = "number")
  xp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x
  xq <- ggplot2::ggplot_build(q)$layout$panel_params[[1]]$x
  expect_equal(xp$get_breaks(), xq$get_breaks())
  ok <- !is.na(xq$get_breaks())
  expect_equal(xp$get_labels()[ok], paste0(100 * as.numeric(xq$get_labels()[ok]), "%"))
})

test_that("visualizeTree() axis argument is partially matched", {
  p <- visualizeTree(Halo_Tree, axis = "perc")
  expect_match(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels(), "%$",
               all = FALSE)
  expect_error(visualizeTree(Halo_Tree, axis = "yes"), "should be one of")
  expect_error(visualizeTree(Halo_Tree, axis = "n"), "should be one of")
  expect_error(visualizeTree(Halo_Tree, axis = TRUE), "must be NULL or a character vector")
})

test_that("visualizeTree() axis can be removed", {
  p <- visualizeTree(Halo_Tree, axis = "none")
  expect_s3_class(p$theme$axis.text.x, "ggplot2::element_blank")
})

test_that("visualizeTree() value labels have solid borders", {
  # The labels should not inherit the dashed linetype of ambiguous branches.
  p <- visualizeTree(halo_chr_tree(), "Scrambling_index", trait = "ChromNumber")
  d <- ggplot2::layer_data(p, layer_index(p, "GeomLabel"))
  expect_true(all(d$linetype %in% c("solid", 1)))
})
