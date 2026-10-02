test_that("makeValueTibble() projects a matrix on internal nodes", {
  tree <- makeValueTibble(Halo_Tree, Halo_PercentDiff, colname = "pd")
  expect_true("pd" %in% names(tree))
  expect_true(all(is.na(tree$pd[tree$isTip])))
  expect_false(anyNA(tree$pd[!tree$isTip]))
  for (node in tree$node[!tree$isTip])
    expect_equal(tree$pd[node], extractValues(node, Halo_Tree, Halo_PercentDiff))
})

test_that("makeValueTibble() matches the values stored in Halo_Tree", {
  tree <- makeValueTibble(Halo_Tree, Halo_PercentDiff, colname = "pd")
  expect_equal(tree$pd, Halo_Tree$Percent_difference)
})

test_that("makeValueTibble() returns a constant for a constant matrix", {
  m <- Halo_PercentDiff
  m[] <- 3
  tree <- makeValueTibble(Halo_Tree, m)
  expect_equal(tree$value[!tree$isTip], rep(3, 5))
})

test_that("makeValueTibble() uses the summary function", {
  tree <- makeValueTibble(Halo_Tree, Halo_PercentDiff, fun = max, colname = "pd")
  expect_equal(tree$pd[7], max(Halo_PercentDiff[childSpecies(Halo_Tree, 7)$left,
                                                childSpecies(Halo_Tree, 7)$right]))
})

test_that("makeValueTibble() returns the tree unchanged when the matrix is NULL", {
  expect_identical(makeValueTibble(Halo_Tree, NULL), Halo_Tree)
})
