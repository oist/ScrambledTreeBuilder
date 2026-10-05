test_that("makeTraitTibble() puts the values on the tips only", {
  tree <- makeTraitTibble(Halo_Tree, halo_taxons(), columns = "ChromNumber")
  tips <- tree$isTip
  expect_equal(tree$ChromNumber[tips], halo_taxons()[tree$label[tips], "ChromNumber"])
  expect_true(all(is.na(tree$ChromNumber[!tips])))
  expect_equal(tree, halo_chr_tree())
})

test_that("makeTraitTibble() matches species by name, not by order", {
  taxons <- halo_taxons()[6:1, ]
  expect_equal(makeTraitTibble(Halo_Tree, taxons), makeTraitTibble(Halo_Tree, halo_taxons()))
})

test_that("makeTraitTibble() adds all columns of the taxon table by default", {
  tree <- makeTraitTibble(Halo_Tree, halo_taxons())
  expect_equal(setdiff(names(tree), names(Halo_Tree)), c("Binomial", "ChromNumber"))
  expect_type(tree$Binomial, "character")
  expect_equal(tree$Binomial[tree$label %in% "Haloferax_volcanii"], "Haloferax volcanii")
})

test_that("makeTraitTibble() requires existing columns", {
  expect_error(makeTraitTibble(Halo_Tree, halo_taxons(), columns = c("ChromNumber", "nope")), "nope")
})

test_that("makeTraitTibble() accepts a named vector", {
  chr <- setNames(halo_taxons()$ChromNumber, rownames(halo_taxons()))
  expect_equal(makeTraitTibble(Halo_Tree, chr, colname = "ChromNumber"), halo_chr_tree())
  expect_true("trait" %in% names(makeTraitTibble(Halo_Tree, chr)))
  expect_error(makeTraitTibble(Halo_Tree, unname(chr)), "named")
})

test_that("makeTraitTibble() keeps the type of the values", {
  taxons <- data.frame(row.names = halo_species,
                       Shape = factor(c("rod", "rod", "rod", "disc", "disc", "rod")),
                       Motile = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))
  tree <- makeTraitTibble(Halo_Tree, taxons)
  expect_s3_class(tree$Shape, "factor")
  expect_equal(levels(tree$Shape), c("disc", "rod"))
  expect_type(tree$Motile, "logical")
})

test_that("makeTraitTibble() warns about tips missing from the table", {
  expect_warning(tree <- makeTraitTibble(Halo_Tree, halo_taxons()[-c(1, 6), ]),
                 "Halobacterium_litoreum.*Salarchaeum_japonicum")
  expect_true(is.na(tree$ChromNumber[tree$label %in% "Halobacterium_litoreum"]))
})

test_that("makeTraitTibble() ignores species absent from the tree", {
  taxons <- rbind(halo_taxons(), data.frame(row.names = "Haloarcula_marismortui",
                                            Binomial = "Haloarcula marismortui", ChromNumber = 1))
  expect_no_warning(tree <- makeTraitTibble(Halo_Tree, taxons))
  expect_equal(tree, makeTraitTibble(Halo_Tree, halo_taxons()))
})

test_that("makeTraitTibble() replaces existing columns", {
  tree <- Halo_Tree
  tree$ChromNumber <- 99
  expect_equal(makeTraitTibble(tree, halo_taxons(), "ChromNumber"), halo_chr_tree())
})

test_that("makeTraitTibble() output can be plotted with a trait", {
  p <- makeTraitTibble(Halo_Tree, halo_taxons()) |> visualizeTree(trait = "ChromNumber")
  expect_equal(p$data$ChromNumber, visualizeTree(halo_chr_tree(), trait = "ChromNumber")$data$ChromNumber)
})
