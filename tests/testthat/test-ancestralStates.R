test_that("ancestralStates() keeps tip values and reconstructs internal nodes", {
  tree <- ancestralStates(Halo_Tree, c(Halobacterium_litoreum  = 3, Halobacterium_noricense = 3,
                                       Halobacterium_salinarum = 3, Salarchaeum_japonicum  = 5,
                                       Haloferax_mediterranei  = 2, Haloferax_volcanii     = 2))
  expect_equal(tree$state[tree$isTip], c(3, 3, 3, 2, 2, 5))
  expect_equal(tree$state[tree$node == 8],  2)  # Haloferax
  expect_equal(tree$state[tree$node == 10], 3)  # Halobacterium
  expect_equal(tree$state[tree$node == 11], 3)
  expect_true(is.na(tree$state[tree$node == 9]))  # 3 or 5
  expect_true(is.na(tree$state[tree$node == 7]))  # root: 2, 3 or 5
})

test_that("ancestralStates() reads values from a taxon table column", {
  tree <- ancestralStates(Halo_Tree, halo_taxons(), column = "ChromNumber")
  expect_true("ChromNumber" %in% names(tree))
  expect_equal(tree$ChromNumber[tree$node == 10], 3)
  expect_equal(tree$ChromNumber[tree$isTip], c(3, 3, 3, 2, 2, 5))
})

test_that("ancestralStates() gives the same state everywhere when all tips agree", {
  tree <- ancestralStates(Halo_Tree, setNames(rep(1, 6), halo_species))
  expect_equal(tree$state, rep(1, 11))
})

test_that("ancestralStates() resolves ambiguous nodes from their parent", {
  # ((A:12, B:12), (C:12, D:14)): the root is 12, so node CD is 12 too and
  # the change to 14 happens on the branch leading to D.
  tree <- ancestralStates(abcd_tree(), c(A = 12, B = 12, C = 12, D = 14))
  expect_equal(tree$state[tree$node == mrca_of(tree, "C", "D")], 12)
  expect_equal(tree$state[tree$node == mrca_of(tree, "A", "D")], 12)
  expect_equal(tree$state[tree$label %in% "D"], 14)
})

test_that("ancestralStates() leaves truly ambiguous nodes as NA", {
  tree <- ancestralStates(abcd_tree(), c(A = 1, B = 2, C = 3, D = 4))
  expect_true(all(is.na(tree$state[!tree$isTip])))
  expect_equal(tree$state[tree$isTip], c(1, 2, 3, 4))
})

test_that("ancestralStates() keeps the type of the values", {
  tree <- ancestralStates(abcd_tree(), c(A = "x", B = "x", C = "y", D = "y"))
  expect_type(tree$state, "character")
  expect_equal(tree$state[tree$node == mrca_of(tree, "C", "D")], "y")
  tree <- ancestralStates(abcd_tree(), c(A = 1L, B = 1L, C = 1L, D = 2L))
  expect_type(tree$state, "integer")
})

test_that("ancestralStates() treats missing tip values as unknown", {
  values <- c(A = 12, B = 12, C = 14)
  expect_warning(tree <- ancestralStates(abcd_tree(), values), "D")
  expect_true(is.na(tree$state[tree$label %in% "D"]))
  expect_equal(tree$state[tree$node == mrca_of(tree, "C", "D")], 14)
  expect_warning(tree <- ancestralStates(abcd_tree(), c(values, D = NA)), "D")
  expect_true(is.na(tree$state[tree$label %in% "D"]))
})

test_that("ancestralStates() ignores values for species absent from the tree", {
  expect_equal(ancestralStates(abcd_tree(), c(A = 1, B = 1, C = 2, D = 2, E = 3)),
               ancestralStates(abcd_tree(), c(A = 1, B = 1, C = 2, D = 2)))
})

test_that("ancestralStates() requires named values or a valid column", {
  expect_error(ancestralStates(abcd_tree(), c(1, 1, 2, 2)), "named")
  expect_error(ancestralStates(Halo_Tree, halo_taxons(), column = "nope"), "nope")
})

test_that("ancestralStates(method = 'ML') reconstructs well-supported clades", {
  tc <- two_clades()
  tree <- ancestralStates(tc$tree, tc$values, method = "ML")
  clade1 <- mrca_of(tree, "sp1", "sp10")
  clade2 <- mrca_of(tree, "sp11", "sp20")
  expect_equal(tree$state[tree$node == clade1], 20)
  expect_equal(tree$state[tree$node == clade2], 22)
  expect_gt(tree$state_prob[tree$node == clade1], 0.95)
  expect_equal(tree$state_prob[tree$isTip], rep(1, 20))
  expect_equal(tree$state[tree$isTip], unname(tc$values[tree$label[tree$isTip]]))
})

test_that("ancestralStates(method = 'ML') leaves poorly supported nodes as NA", {
  tc <- two_clades()
  tree <- ancestralStates(tc$tree, tc$values, method = "ML")
  root <- tree$node[tree$parent == tree$node]
  expect_true(is.na(tree$state[tree$node == root]))
  expect_equal(tree$state_prob[tree$node == root], 0.5, tolerance = 1e-3)
})

test_that("ancestralStates(method = 'ML') uses the threshold", {
  tree <- ancestralStates(Halo_Tree, halo_taxons(), "ChromNumber", method = "ML")
  expect_true(all(is.na(tree$ChromNumber[!tree$isTip])))
  expect_true(all(tree$ChromNumber_prob[!tree$isTip] < 0.95))
  tree <- ancestralStates(Halo_Tree, halo_taxons(), "ChromNumber", method = "ML", threshold = 0.8)
  expect_equal(tree$ChromNumber[tree$node == 10], 3)
  expect_equal(tree$ChromNumber[tree$node == 11], 3)
  expect_true(is.na(tree$ChromNumber[tree$node == 8]))  # 2, with probability 0.77
  expect_equal(tree$ChromNumber_prob[tree$node == 8], 0.77, tolerance = 0.01)
})

test_that("ancestralStates(model = 'ordered') allows unobserved intermediate states", {
  tc <- two_clades()
  tree <- ancestralStates(tc$tree, tc$values, method = "ML", model = "ordered", threshold = 0.4)
  root <- tree$node[tree$parent == tree$node]
  expect_equal(tree$state[tree$node == root], 21)
  tree <- ancestralStates(tc$tree, tc$values, method = "ML", model = "ER", threshold = 0.4)
  expect_true(tree$state[tree$node == root] %in% c(20, 22))
})

test_that("ancestralStates(model = 'ordered') requires integer values", {
  expect_error(ancestralStates(abcd_tree(), c(A = "x", B = "x", C = "y", D = "y"),
                               method = "ML", model = "ordered"), "integer")
})

test_that("ancestralStates(method = 'ML') warns when the model can not be fitted", {
  # Values alternating within sister pairs carry no phylogenetic signal.
  expect_warning(ancestralStates(abcd_tree(), c(A = 1, B = 3, C = 1, D = 3),
                                 method = "ML", model = "ordered"),
                 "uniform")
})

test_that("ancestralStates(method = 'ML') needs a value for every tip", {
  expect_error(ancestralStates(abcd_tree(), c(A = 12, B = 12, C = 14), method = "ML"), "D")
})

test_that("ancestralStates() rejects unknown methods", {
  expect_error(ancestralStates(abcd_tree(), c(A = 1, B = 1, C = 1, D = 1), method = "magic"))
})

test_that("ancestralStates(method = 'ML') does not depend on the scale of branch lengths", {
  # Regression: depending on the scale of the branch lengths, ape::ace()
  # sometimes failed to fit the ordered model ("NA/NaN/Inf in foreign
  # function call"), for instance here with 20 vs 25 and a 10000-fold
  # shrink, or 20 vs 30 and a 100-fold shrink.
  tc <- two_clades()
  for (high in c(25, 30)) {
    values <- setNames(ifelse(tc$values == 20, 20, high), names(tc$values))
    fits <- lapply(c(1, 1e2, 1e4), \(scale) {
      tree <- tc$tree
      tree$branch.length <- tree$branch.length / scale
      ancestralStates(tree, values, method = "ML", model = "ordered", threshold = 0.4)
    })
    expect_equal(fits[[2]][, c("state", "state_prob")], fits[[1]][, c("state", "state_prob")], tolerance = 1e-4)
    expect_equal(fits[[3]][, c("state", "state_prob")], fits[[1]][, c("state", "state_prob")], tolerance = 1e-4)
  }
})
