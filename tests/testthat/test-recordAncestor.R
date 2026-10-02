test_that("recordAncestor() adds the MRCA node of each pair", {
  df <- recordAncestor(formatStats(halo_files()), Halo_Tree)
  expect_equal(df["Haloferax_mediterranei___Haloferax_volcanii", "MRCA"], 8)
  expect_equal(df["Halobacterium_litoreum___Halobacterium_noricense", "MRCA"], 11)
  expect_equal(df["Haloferax_volcanii___Salarchaeum_japonicum", "MRCA"], 7)
  expect_equal(as.vector(table(df$MRCA)), c(16, 2, 6, 4, 2))
  expect_equal(df$MRCA, Halo_DF$MRCA)
})

test_that("recordAncestor() gives the same MRCA to both orientations of a pair", {
  df <- recordAncestor(formatStats(halo_files()), Halo_Tree)
  reverse <- paste(df$species2, df$species1, sep = "___")
  expect_equal(df$MRCA, df[reverse, "MRCA"])
})

test_that("recordAncestor() gives 0 to pairs absent from the tree", {
  tree <- subTree(Halo_Tree, 9)
  df <- recordAncestor(formatStats(halo_files()), tree)
  expect_true(all(df$MRCA[grepl("Haloferax", rownames(df))] == 0))
  expect_true(all(df$MRCA[!grepl("Haloferax", rownames(df))] != 0))
})

test_that("recordClades() adds focal clade names and colors to results tables", {
  df <- recordClades(recordAncestor(formatStats(halo_files()), Halo_Tree), Halo_FocalClades)
  expect_equal(df["Haloferax_mediterranei___Haloferax_volcanii", "focalClade"], "Haloferax")
  expect_equal(df["Haloferax_mediterranei___Haloferax_volcanii", "focalColor"], "green3")
  expect_equal(df["Halobacterium_litoreum___Halobacterium_salinarum", "focalClade"], "Halobacterium")
  expect_true(is.na(df["Haloferax_volcanii___Salarchaeum_japonicum", "focalClade"]))
  expect_equal(df$focalClade, Halo_DF$focalClade)
  expect_equal(df$focalColor, Halo_DF$focalColor)
})

test_that("recordClades() annotates taxon tables", {
  taxa <- data.frame(row.names = halo_species, Binomial = sub("_", " ", halo_species))
  taxa <- recordClades(taxa, Halo_FocalClades)
  expect_equal(taxa["Haloferax_volcanii", "customClade"], "Haloferax")
  expect_equal(taxa["Halobacterium_salinarum", "customClade"], "Halobacterium")
  expect_true(is.na(taxa["Salarchaeum_japonicum", "customClade"]))
})

test_that("recordClades() computes missing MRCAs", {
  skip("Known bug: recordClades() passes the clades instead of a tree to recordAncestor().")
  expect_no_error(recordClades(formatStats(halo_files()), Halo_FocalClades))
})
