test_that("eog() builds a command with brace expansions for a clade", {
  expect_equal(eog(Halo_Tree, 8),
               "eog results_{Haloferax_mediterranei,Haloferax_volcanii}/alignment/*___{Haloferax_mediterranei,Haloferax_volcanii}.o2o_plt.png")
  expect_equal(eog(Halo_Tree, Halo_FocalClades$Haloferax), eog(Halo_Tree, 8))
})

test_that("eog() can restrict to target genomes", {
  expect_equal(eog(Halo_Tree, 8, target = "volc"),
               "eog results_Haloferax_volcanii/alignment/*___{Haloferax_mediterranei,Haloferax_volcanii}.o2o_plt.png")
  expect_match(eog(Halo_Tree, 10, target = "Halobacterium_[ln]"),
               "^eog results_\\{Halobacterium_litoreum,Halobacterium_noricense\\}/")
  expect_error(eog(Halo_Tree, 8, target = "no_such_species"), "does not match")
})
