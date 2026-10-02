test_that("subTree() extracts the clade under a node", {
  s <- subTree(Halo_Tree, 9)
  expect_s3_class(s, "tbl_tree")
  expect_setequal(s$label[s$isTip], c("Halobacterium_litoreum", "Halobacterium_noricense",
                                      "Halobacterium_salinarum", "Salarchaeum_japonicum"))
  expect_equal(nrow(s), 7)
  expect_equal(sort(s$y[s$isTip]), 1:4)
})

test_that("subTree() records the original node IDs", {
  s <- subTree(Halo_Tree, 9)
  expect_setequal(s$node.orig, c(1, 2, 3, 6, 9, 10, 11))
  tip <- s$label %in% "Salarchaeum_japonicum"
  expect_equal(s$node.orig[tip], Halo_Tree$node[Halo_Tree$label %in% "Salarchaeum_japonicum"])
})

test_that("subTree() accepts focal clades", {
  expect_equal(subTree(Halo_Tree, Halo_FocalClades$Haloferax), subTree(Halo_Tree, 8))
})
