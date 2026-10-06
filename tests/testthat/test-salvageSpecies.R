# Two caterpillar clades A and B, with nodes at distances 1, 2 and 3, joined
# at distance 10.  Terminal branches are half the distance to the closest
# relative: 0.5 for A1, A2, 1 for A3 and 1.5 for A4 (same for B).
twin <- function() {
  sp <- c(paste0("A", 1:4), paste0("B", 1:4))
  cat4 <- matrix(c(0, 1, 2, 3,
                   1, 0, 2, 3,
                   2, 2, 0, 3,
                   3, 3, 3, 0), 4)
  D <- matrix(10, 8, 8, dimnames = list(sp, sp))
  D[1:4, 1:4] <- cat4
  D[5:8, 5:8] <- cat4
  tree <- makeTidyTree(D) |> makeValueTibble(D, colname = "d")
  tree$bootstrap <- ifelse(tree$isTip, NA, 100)
  tree
}
node_at <- function(tree, a, b) mrca_of(tree, a, b)

test_that("salvageSpecies() returns the species on the longest terminal branches", {
  tree <- twin()
  s <- salvageSpecies(tree, "d", long = 2, support = NULL)
  expect_setequal(s$species, c("A4", "B4"))
  expect_match(s$reason, "long terminal branch")
  bl <- setNames(tree$branch.length, tree$label)
  expect_true(all(bl[s$species] == max(bl[tree$isTip])))
})

test_that("salvageSpecies() keeps both sides of weakly supported deep nodes", {
  tree <- twin()
  root <- node_at(tree, "A1", "B1")
  a3   <- node_at(tree, "A1", "A4")
  tree$bootstrap[tree$node %in% c(root, a3)] <- 30
  s <- salvageSpecies(tree, "d", long = 0, support = 70, depth = 2.5)
  # Deepest first: the root needs one A and one B (the longest branches, A4
  # and B4); then node A1-A4 already has A4 and needs one of A1-A3 (A3).
  expect_setequal(s$species[1:2], c("A4", "B4"))
  expect_equal(s$species[3], "A3")
  expect_match(s$reason, "weakly supported node")
  for (n in c(root, a3)) {
    sides <- childSpecies(tree, n)
    expect_true(any(sides$left %in% s$species) && any(sides$right %in% s$species))
  }
})

test_that("salvageSpecies() ignores well supported and shallow nodes", {
  tree <- twin()
  tree$bootstrap[tree$node %in% node_at(tree, "A1", "A2")] <- 10  # at distance 1
  expect_equal(nrow(salvageSpecies(tree, "d", long = 0, support = 70, depth = 2.5)), 0)
  expect_equal(nrow(salvageSpecies(tree, "d", long = 0, support = 5, depth = 0)), 0)
  expect_equal(salvageSpecies(tree, "d", long = 0, support = 70, depth = 0)$species, c("A1", "A2"))
})

test_that("salvageSpecies() uses the species already kept", {
  tree <- twin()
  tree$bootstrap[tree$node %in% node_at(tree, "A1", "B1")] <- 30
  s <- salvageSpecies(tree, "d", long = 0, support = 70, depth = 5, keep = c("A1", "B2"))
  expect_equal(nrow(s), 0)
  s <- salvageSpecies(tree, "d", long = 0, support = 70, depth = 5, keep = "A1")
  expect_equal(s$species, "B4")
})

test_that("salvageSpecies() deepens by default to half the largest node value", {
  tree <- twin()
  tree$bootstrap[tree$node %in% node_at(tree, "A1", "A4")] <- 30  # at distance 3
  expect_equal(nrow(salvageSpecies(tree, "d", long = 0)), 0)          # 3 < 10 / 2
  tree$bootstrap[tree$node %in% node_at(tree, "A1", "B1")] <- 30  # at distance 10
  expect_setequal(salvageSpecies(tree, "d", long = 0)$species, c("A4", "B4"))
})

test_that("salvageSpecies() warns when the tree has no bootstrap support", {
  expect_warning(s <- salvageSpecies(Halo_Tree, "Percent_difference", long = 1), "bootstrap")
  expect_equal(nrow(s), 1)
  expect_equal(s$species, Halo_Tree$label[which.max(ifelse(Halo_Tree$isTip, Halo_Tree$branch.length, NA))])
  expect_no_warning(salvageSpecies(Halo_Tree, "Percent_difference", long = 1, support = NULL))
})

test_that("salvageSpecies() output can be kept by subsampleSpecies()", {
  tree <- twin()
  tree$bootstrap[tree$node %in% node_at(tree, "A1", "A4")] <- 30
  s <- salvageSpecies(tree, "d", long = 1, support = 70, depth = 2.5)
  set.seed(1)
  sub <- subsampleSpecies(tree, "d", n = 4, keep = s$species)
  expect_true(all(s$species %in% sub@species))
})
