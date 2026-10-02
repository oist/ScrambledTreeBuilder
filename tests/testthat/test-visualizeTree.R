test_that("visualizeTree() labels internal nodes with their IDs by default", {
  p <- visualizeTree(Halo_Tree)
  expect_s3_class(p, "ggtree")
  b <- ggplot2::ggplot_build(p)
  labels <- b$data[[length(b$data)]]$label
  expect_setequal(labels, Halo_Tree$node)
  expect_equal(p$theme$legend.position, "none")
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
