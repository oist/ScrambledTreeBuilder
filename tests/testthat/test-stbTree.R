test_that("orderWithTree() sorts the matrix like the plotted tree leaves", {
  m <- orderWithTree(Halo_PercentDiff, Halo_Tree)
  expect_equal(rownames(m), colnames(m))
  tips <- Halo_Tree[Halo_Tree$isTip, ]
  expect_equal(rownames(m), tips$label[order(tips$y, decreasing = TRUE)])
  expect_equal(m["Haloferax_volcanii", "Salarchaeum_japonicum"],
               Halo_PercentDiff["Haloferax_volcanii", "Salarchaeum_japonicum"])
})
