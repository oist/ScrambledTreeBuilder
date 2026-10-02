test_that("extractValues() summarises the left x right block of a node", {
  m <- Halo_PercentDiff
  expect_equal(extractValues(8, Halo_Tree, m),
               m["Haloferax_mediterranei", "Haloferax_volcanii"])
  ch <- childSpecies(Halo_Tree, 7)
  expect_equal(extractValues(7, Halo_Tree, m), mean(m[ch$left, ch$right]))
  expect_equal(extractValues(7, Halo_Tree, m, fun = max), max(m[ch$left, ch$right]))
})
