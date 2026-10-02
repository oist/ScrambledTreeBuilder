# Consistency of the lazy-loaded example data.  If these fail, re-run the
# scripts in data-raw/.

test_that("Halo_* objects are consistent with each other", {
  expect_equal(nrow(Halo_DF), 30)
  expect_equal(rownames(Halo_PercentDiff), halo_species)
  expect_setequal(Halo_Tree$label[Halo_Tree$isTip], halo_species)
  for (cl in Halo_FocalClades)
    expect_true(all(cl@genomeIDs %in% halo_species))
})

test_that("oikData tables can be overlaid on MRCA_2D_plot() data", {
  expect_named(oikData, c("2025_02_17", "2025_02_25", "2025_07_17"))
  plot_data <- MRCA_2D_plot(Halo_DF, Halo_FocalClades)$data
  expect_equal(names(oikData[["2025_07_17"]]), names(plot_data))
  expect_no_error(rbind(plot_data, oikData[["2025_07_17"]]))
  for (tb in oikData)
    expect_true(all(c("MRCA", "x", "y", "xerr", "yerr", "n", "clade") %in% names(tb)))
})
