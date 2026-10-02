test_that("subMatrix() subsets a matrix to the species of a clade", {
  m <- subMatrix(Halo_Tree, Halo_PercentDiff, 10)
  sp <- c("Halobacterium_litoreum", "Halobacterium_noricense", "Halobacterium_salinarum")
  expect_setequal(rownames(m), sp)
  expect_equal(m, Halo_PercentDiff[rownames(m), rownames(m)])
})

test_that("subMatrix() accepts focal clades and subtrees", {
  expect_equal(subMatrix(Halo_Tree, Halo_PercentDiff, Halo_FocalClades$Halobacterium),
               subMatrix(Halo_Tree, Halo_PercentDiff, 10))
  expect_equal(subTree(Halo_Tree, 10) |> subMatrix(Halo_PercentDiff),
               subMatrix(Halo_Tree, Halo_PercentDiff, 10))
})

test_that("subMatrix() without a clade keeps all the species", {
  expect_setequal(rownames(subMatrix(Halo_Tree, Halo_PercentDiff)), halo_species)
})

test_that("subMatrix() can abbreviate genus names", {
  m <- subMatrix(Halo_Tree, Halo_PercentDiff, simplify.names = TRUE)
  expect_true("H_volcanii"  %in% rownames(m))
  expect_true("S_japonicum" %in% colnames(m))
  expect_equal(rownames(m), colnames(m))
})
