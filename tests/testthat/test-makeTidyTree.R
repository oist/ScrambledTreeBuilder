test_that("makeTidyTree() builds a stbTree with tips and internal nodes", {
  tree <- makeTidyTree(halo_sym())
  expect_s3_class(tree, "tbl_tree")
  expect_equal(nrow(tree), 11)
  expect_equal(sum(tree$isTip), 6)
  expect_setequal(tree$label[tree$isTip], halo_species)
  expect_true(all(is.na(tree$label[!tree$isTip])))
  expect_true(all(is.na(tree$bootstrap)))
  expect_true(all(c("parent", "node", "branch.length", "label", "isTip", "y", "bootstrap") %in% names(tree)))
})

test_that("makeTidyTree() reproduces the topology of Halo_Tree", {
  tree <- makeTidyTree(halo_sym())
  expect_equal(tree[, c("parent", "node", "label", "isTip", "y")],
               Halo_Tree[, c("parent", "node", "label", "isTip", "y")])
  expect_equal(tree$branch.length, Halo_Tree$branch.length)
})

test_that("makeTidyTree() clusters with UPGMA", {
  tree <- makeTidyTree(halo_sym())
  # The closest pair (Haloferax) joins with half its distance as branch length.
  d <- halo_sym()["Haloferax_mediterranei", "Haloferax_volcanii"]
  expect_equal(tree$branch.length[tree$label %in% "Haloferax_volcanii"], d / 2)
})

test_that("makeTidyTree() warns on asymmetric matrices", {
  expect_warning(makeTidyTree(Halo_PercentDiff), "not symmetric")
})

test_that("makeTidyTree() records bootstrap support on internal nodes", {
  set.seed(1)
  tree <- makeTidyTree(halo_sym(), n_bootstrap = 20)
  expect_true(all(is.na(tree$bootstrap[tree$isTip])))
  expect_true(all(tree$bootstrap[!tree$isTip] >= 0 & tree$bootstrap[!tree$isTip] <= 100))
  expect_true(all(is.na(tree$label[!tree$isTip])))
})

test_that("makeItConvenient() recomputes isTip and y", {
  tree <- Halo_Tree
  tree$isTip <- NULL
  tree$y     <- NULL
  tree <- ScrambledTreeBuilder:::makeItConvenient(tree)
  expect_equal(tree$isTip, Halo_Tree$isTip)
  expect_equal(tree$y,     Halo_Tree$y)
})
