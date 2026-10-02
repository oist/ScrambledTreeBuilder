test_that("focalClade() records the clade at the MRCA of two species", {
  cl <- focalClade(Halo_Tree, "Halobacterium_noricense", "Halobacterium_salinarum",
                   "blue", "Halobacterium")
  expect_true(S7::S7_inherits(cl, ScrambledTreeBuilder:::FocalClade))
  expect_equal(cl@nodeID, 10)
  expect_setequal(cl@genomeIDs, c("Halobacterium_litoreum", "Halobacterium_noricense",
                                  "Halobacterium_salinarum"))
  expect_setequal(cl@nodeList, c(1, 2, 3, 10, 11))
  expect_equal(cl@color, "blue")
  expect_equal(cl@displayName, "Halobacterium")
})

test_that("focalClade() objects match the Halo_FocalClades example data", {
  cl <- focalClade(Halo_Tree, "Haloferax_mediterranei", "Haloferax_volcanii",
                   "green3", "Haloferax")
  expect_equal(cl, Halo_FocalClades$Haloferax)
})

test_that("FocalCladeList() is a named list of focal clades", {
  expect_true(S7::S7_inherits(Halo_FocalClades, FocalCladeList))
  expect_named(Halo_FocalClades, c("Halobacterium", "Haloferax"))
  expect_equal(Halo_FocalClades$Haloferax@nodeID, 8)
})

test_that("Focal clades print one line per clade", {
  expect_equal(capture.output(print(Halo_FocalClades)),
               c("Halobacterium, node ID: 10, number of genomes: 3",
                 "Haloferax, node ID: 8, number of genomes: 2"))
})

test_that("Adding focal clades to a tree plot adds highlight layers", {
  p <- visualizeTree(Halo_Tree)
  expect_length((p + Halo_FocalClades$Haloferax)$layers, length(p$layers) + 1)
  expect_length((p + Halo_FocalClades)$layers,           length(p$layers) + 2)
  expect_s3_class(p + Halo_FocalClades, "ggtree")
})
