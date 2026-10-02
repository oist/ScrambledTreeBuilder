test_that("treeHeatMap() returns a pheatmap with or without tree and clades", {
  expect_s3_class(treeHeatMap(Halo_PercentDiff, silent = TRUE), "pheatmap")
  expect_s3_class(treeHeatMap(Halo_PercentDiff, Halo_Tree, silent = TRUE), "pheatmap")
  expect_s3_class(treeHeatMap(Halo_PercentDiff, Halo_Tree, Halo_FocalClades, silent = TRUE), "pheatmap")
})

test_that("treeHeatMap() sorts rows like the tree", {
  h <- treeHeatMap(Halo_PercentDiff, Halo_Tree, silent = TRUE)
  row_names <- h$gtable$grobs[[which(h$gtable$layout$name == "row_names")]]$label
  expect_equal(row_names, rownames(orderWithTree(Halo_PercentDiff, Halo_Tree)))
})

test_that("treeHeatMap() draws a clade annotation", {
  h <- treeHeatMap(Halo_PercentDiff, Halo_Tree, Halo_FocalClades, silent = TRUE)
  expect_true("row_annotation" %in% h$gtable$layout$name)
})
