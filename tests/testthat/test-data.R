# Consistency of the lazy-loaded example data.  If these fail, re-run the
# scripts in data-raw/.

test_that("Halo_* objects are consistent with each other", {
  expect_equal(nrow(Halo_DF), 30)
  expect_equal(rownames(Halo_PercentDiff), halo_species)
  expect_setequal(Halo_Tree$label[Halo_Tree$isTip], halo_species)
  for (cl in Halo_FocalClades)
    expect_true(all(cl@genomeIDs %in% halo_species))
  expect_equal(rownames(Halo_Taxons), halo_species)
})

test_that("Halo_Taxons illustrates the features of trait plots", {
  tree <- makeTraitTibble(Halo_Tree, Halo_Taxons) |> ancestralStates(Halo_Taxons, "Toy_trait")
  state <- \(a, b) tree$Toy_trait[tree$node == mrca_of(tree, a, b)]
  # One value per focal clade, with an exception on a terminal branch.
  expect_equal(state("Halobacterium_litoreum", "Halobacterium_salinarum"), 3)
  expect_equal(state("Halobacterium_litoreum", "Halobacterium_noricense"), 3)
  expect_equal(Halo_Taxons["Halobacterium_litoreum", "Toy_trait"], 4)
  expect_equal(state("Haloferax_mediterranei", "Haloferax_volcanii"), 2)
  # Ambiguous ancestors.
  expect_true(is.na(state("Salarchaeum_japonicum", "Halobacterium_salinarum")))
  expect_true(is.na(state("Salarchaeum_japonicum", "Haloferax_volcanii")))
})

test_that("oikData tables can be overlaid on MRCA_2D_plot() data", {
  expect_named(oikData, c("2025_02_17", "2025_02_25", "2025_07_17"))
  plot_data <- MRCA_2D_plot(Halo_DF, Halo_FocalClades)$data
  expect_equal(names(oikData[["2025_07_17"]]), names(plot_data))
  expect_no_error(rbind(plot_data, oikData[["2025_07_17"]]))
  for (tb in oikData)
    expect_true(all(c("MRCA", "x", "y", "xerr", "yerr", "n", "clade") %in% names(tb)))
})
