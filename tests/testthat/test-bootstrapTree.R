test_that("residualBootstrapTree() returns a reference tree and replicates", {
  set.seed(1)
  bt <- residualBootstrapTree(Halo_PercentDiff, n_bootstrap = 20)
  expect_s3_class(bt, "residualBootstrapTreeResult")
  expect_named(bt, c("tree", "bootstrap_trees", "support", "call"))
  expect_s3_class(bt$tree, "phylo")
  expect_length(bt$bootstrap_trees, 20)
  expect_setequal(bt$tree$tip.label, halo_species)
  expect_length(bt$support, bt$tree$Nnode)
  expect_true(all(bt$support >= 0 & bt$support <= 100))
  expect_equal(bt$tree$node.label, as.character(round(bt$support)))
})

test_that("residualBootstrapTree() is reproducible with a seed", {
  set.seed(2); a <- residualBootstrapTree(Halo_PercentDiff, n_bootstrap = 10)
  set.seed(2); b <- residualBootstrapTree(Halo_PercentDiff, n_bootstrap = 10)
  expect_equal(a$support, b$support)
})

test_that("residualBootstrapTree() validates its input", {
  m <- Halo_PercentDiff
  expect_error(residualBootstrapTree(unname(m), 2), "names")
  expect_error(residualBootstrapTree(m[, rev(colnames(m))], 2), "identical")
  m[1, 2] <- NA
  expect_error(residualBootstrapTree(m, 2), "NA")
})
