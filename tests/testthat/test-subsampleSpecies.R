# A random ultrametric tree of 40 species, with many shallow nodes, and the
# pairwise distances stored on its internal nodes.
coal_tree <- function(n = 40) {
  set.seed(1)
  D <- ape::cophenetic.phylo(ape::rcoal(n, tip.label = paste0("sp", seq_len(n))))
  makeTidyTree(D) |> makeValueTibble(D, colname = "d")
}

# Two clades with the same shape: caterpillars with nodes at distances 1, 2
# and 3, joined at distance 10.
twin_tree <- function() {
  sp <- c(paste0("A", 1:4), paste0("B", 1:4))
  cat4 <- matrix(c(0, 1, 2, 3,
                   1, 0, 2, 3,
                   2, 2, 0, 3,
                   3, 3, 3, 0), 4)
  D <- matrix(10, 8, 8, dimnames = list(sp, sp))
  D[1:4, 1:4] <- cat4
  D[5:8, 5:8] <- cat4
  makeTidyTree(D) |> makeValueTibble(D, colname = "d")
}
twin_clades <- function(tree) {
  FocalCladeList(A = focalClade(tree, "A1", "A4", "red",  "A"),
                 B = focalClade(tree, "B1", "B4", "blue", "B"))
}

# MRCAs of all pairs of species, computed independently.
pair_mrcas <- function(tree, species) {
  pairs <- utils::combn(species, 2)
  sort(unique(apply(pairs, 2, \(p) mrca_of(tree, p[1], p[2]))))
}

test_that("subsampleSpecies() keeps n species and their n - 1 MRCAs", {
  tree <- coal_tree()
  set.seed(1)
  s <- subsampleSpecies(tree, "d", n = 12)
  expect_length(s@species, 12)
  expect_true(all(s@species %in% tree$label[tree$isTip]))
  expect_equal(sort(s@nodes), pair_mrcas(tree, s@species))
})

test_that("subsampleSpecies() reports the nodes per bin", {
  tree <- coal_tree()
  set.seed(1)
  s <- subsampleSpecies(tree, "d", n = 12, breaks = 4)
  d <- tree$d[!tree$isTip]
  bins <- cut(d, seq(0, max(d), length.out = 5), include.lowest = TRUE)
  expect_equal(s@report$available, as.vector(table(bins)))
  expect_equal(s@report$kept,
               as.vector(table(bins[match(s@nodes, tree$node[!tree$isTip])])))
  expect_equal(sum(s@report$target), 11)
  expect_true(all(s@report$target <= s@report$available))
})

test_that("subsampleSpecies() spreads the MRCAs over the bins", {
  tree <- coal_tree()
  d <- tree$d[!tree$isTip]
  # Most nodes are shallow: a random subsample would keep mostly shallow ones.
  expect_gt(mean(d < max(d) / 4), 0.5)
  set.seed(1)
  s <- subsampleSpecies(tree, "d", n = 10, breaks = 4)
  expect_equal(s@report$kept, s@report$target)
})

test_that("subsampleSpecies() accepts explicit breaks", {
  tree <- twin_tree()
  s <- subsampleSpecies(tree, "d", n = 3, breaks = c(0, 5, 10))
  expect_equal(s@report$available, c(6, 1))
  expect_equal(s@report$kept, c(1, 1))
})

test_that("subsampleSpecies() prefers MRCAs from different clades in a bin", {
  tree <- twin_tree()
  for (seed in 1:5) {
    set.seed(seed)
    s <- subsampleSpecies(tree, "d", n = 4, breaks = c(0, 5, 10), clades = twin_clades(tree))
    shallow <- setdiff(s@nodes, mrca_of(tree, "A1", "B1"))
    sides <- sapply(shallow, \(n) substr(childSpecies(tree, n)$left[1], 1, 1))
    expect_setequal(sides, c("A", "B"))
  }
})

test_that("subsampleSpecies() finds n from the number of replicates per bin", {
  tree <- coal_tree()
  set.seed(1)
  s <- subsampleSpecies(tree, "d", replicates = 3, breaks = 4)
  expect_equal(s@report$target, pmin(s@report$available, 3))
  expect_true(all(s@report$kept >= s@report$target))
  expect_gte(length(s@species), sum(s@report$target) + 1)
  expect_length(s@nodes, length(s@species) - 1)
})

test_that("waterFill() spreads nodes evenly and keeps the span", {
  available <- c(10, 1, 10, 2, 1)
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 9),  c(2, 1, 3, 2, 1))
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 24), available)
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 1),  c(0, 0, 0, 0, 1))
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 2),  c(1, 0, 0, 0, 1))
})

test_that("waterFill() starts from a floor", {
  available <- c(10, 1, 10, 2, 1)
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 9, floor = c(0, 0, 0, 0, 1)), c(2, 1, 3, 2, 1))
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 9, floor = c(0, 0, 0, 2, 1)), c(2, 1, 3, 2, 1))
  expect_equal(ScrambledTreeBuilder:::waterFill(available, 3, floor = c(0, 0, 0, 2, 1)), c(0, 0, 0, 2, 1))
})

test_that("subsampleSpecies() reaches empty lineages behind an overfilled bin", {
  # Clades A and B have nodes at 1, 2 and 3, clade W at 2, 3 and 4.  Deep
  # nodes: X-Y at 6, W at 8, Z at 9, A-B at 10 and the root at 12.  The
  # species kept fill the deep bin over target, and the shallowest bin can
  # only be completed with two species of B, whose join node is deep.
  # Species of W add nodes to bins already at target.
  cat4 <- \(h) { m <- matrix(h[3], 4, 4); m[1:3, 1:3] <- h[2]; m[1:2, 1:2] <- h[1]; diag(m) <- 0; m }
  sp <- c(paste0("A", 1:4), paste0("B", 1:4), paste0("W", 1:4), "X", "Y", "Z")
  D <- matrix(12, 15, 15, dimnames = list(sp, sp))
  D[1:8, 1:8] <- 10
  D[1:4, 1:4] <- cat4(1:3)
  D[5:8, 5:8] <- cat4(1:3)
  D[9:15, 9:15] <- 9
  D[9:14, 9:14] <- 8
  D[9:12, 9:12] <- cat4(2:4)
  D[13:14, 13:14] <- matrix(c(0, 6, 6, 0), 2)
  diag(D) <- 0
  tree <- makeTidyTree(D) |> makeValueTibble(D, colname = "d")
  for (seed in 1:5) {
    set.seed(seed)
    s <- subsampleSpecies(tree, "d", replicates = 2, breaks = c(0, 1.5, 5, 12),
                          keep = c(paste0("A", 1:4), "W1", "X", "Y", "Z"))
    expect_equal(s@report$kept[1], 2)
    expect_length(s@species, 10)
  }
})

test_that("subsampleSpecies() counts the nodes of kept species in the targets", {
  tree <- coal_tree()
  set.seed(2)
  keep <- sample(tree$label[tree$isTip], 6)
  for (seed in 1:3) {
    set.seed(seed)
    s <- subsampleSpecies(tree, "d", n = 14, breaks = 4, keep = keep)
    expect_equal(sum(pmax(s@report$target - s@report$kept, 0)), 0)
    expect_true(all(keep %in% s@species))
  }
})

test_that("subsampleSpecies() keeps them when only one species is free", {
  # sample(x) on a single number x permutes 1:x.
  tree <- coal_tree()
  for (seed in 1:5) {
    set.seed(seed)
    keep <- sample(tree$label[tree$isTip], 5)
    expect_true(all(keep %in% subsampleSpecies(tree, "d", n = 6, keep = keep)@species))
  }
})

test_that("subsampleSpecies() ignores duplicates in the species to keep", {
  tree <- coal_tree()
  set.seed(1)
  s <- subsampleSpecies(tree, "d", n = 8, keep = c("sp1", "sp2", "sp1"))
  expect_length(s@species, 8)
  expect_false(anyDuplicated(s@species) > 0)
  expect_length(s@nodes, 7)
})

test_that("subsampleSpecies() keeps the species it is asked to keep", {
  tree <- coal_tree()
  set.seed(1)
  s <- subsampleSpecies(tree, "d", n = 8, keep = c("sp1", "sp2", "sp3"))
  expect_true(all(c("sp1", "sp2", "sp3") %in% s@species))
  expect_length(s@species, 8)
})

test_that("subsampleSpecies() handles the extreme sizes", {
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 6)
  expect_setequal(s@species, Halo_Tree$label[Halo_Tree$isTip])
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 2)
  root <- Halo_Tree$node[Halo_Tree$parent == Halo_Tree$node]
  expect_equal(s@nodes, root)
})

test_that("subsampleSpecies() accepts a vector of node values", {
  set.seed(1)
  a <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4)
  set.seed(1)
  b <- subsampleSpecies(Halo_Tree, Halo_Tree$Percent_difference, n = 4)
  expect_equal(a@species, b@species)
})

test_that("subsampleSpecies() checks its arguments", {
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference"), "n.*replicates")
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference", n = 3, replicates = 2), "n.*replicates")
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference", n = 1), "between 2 and 6")
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference", n = 7), "between 2 and 6")
  expect_error(subsampleSpecies(Halo_Tree, "nope", n = 3), "nope")
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference", n = 3, keep = "nope"), "nope")
  expect_error(subsampleSpecies(Halo_Tree, "Percent_difference", n = 2,
                                keep = halo_species[1:3]), "keep")
})

test_that("SpeciesSubsample objects print a summary", {
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4)
  expect_output(print(s), "4 of 6 species")
})

test_that("plotSubsample() draws the kept species and MRCAs on the tree", {
  set.seed(1)
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4, clades = Halo_FocalClades)
  p <- plotSubsample(s, panel = "tree")
  expect_s3_class(p, "ggtree")
  i <- which(sapply(p$layers, \(l) identical(l$aes_params$fill, "gold")))
  expect_length(i, 1)
  expect_setequal(ggplot2::layer_data(p, i)$y, p$data$y[p$data$node %in% s@nodes])
})

test_that("plotSubsample() counts the MRCAs per bin", {
  set.seed(1)
  s <- subsampleSpecies(coal_tree(), "d", n = 10, breaks = 4)
  p <- plotSubsample(s, panel = "bins")
  d <- ggplot2::layer_data(p, 1)
  counts <- as.vector(tapply(d$count, d$PANEL, sum))
  expect_equal(counts, c(39, 9))
})

test_that("plotSubsample() shows which pairs are kept", {
  set.seed(1)
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4)
  p <- plotSubsample(s, Halo_DF, x = "percent_difference_global", panel = "pairs")
  kept <- Halo_DF$species1 %in% s@species & Halo_DF$species2 %in% s@species
  expect_equal(nrow(ggplot2::layer_data(p, 2)), sum(kept))
  expect_equal(nrow(ggplot2::layer_data(p, 1)), sum(!kept))
  expect_error(plotSubsample(s, panel = "pairs"), "pairwise_data")
})

test_that("plotSubsample() combines the panels", {
  skip_if_not_installed("patchwork")
  s <- subsampleSpecies(Halo_Tree, "Percent_difference", n = 4)
  expect_s3_class(plotSubsample(s), "patchwork")
  expect_s3_class(plotSubsample(s, Halo_DF), "patchwork")
})
