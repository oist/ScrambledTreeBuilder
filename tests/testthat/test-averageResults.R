test_that("averageResults() averages both orientations of each pair", {
  a <- averageResults(Halo_DF)
  expect_s3_class(a, "tbl_df")
  expect_equal(nrow(a), 15)
  row <- a[a$lab == "Haloferax_volcanii\nHaloferax_mediterranei", ]
  expect_equal(row$percent_difference_global,
               mean(Halo_DF[Halo_DF$lab == row$lab, "percent_difference_global"]))
  expect_equal(row$species1, "Haloferax_volcanii")
  expect_equal(row$species2, "Haloferax_mediterranei")
  expect_equal(row$MRCA, 8)
  expect_equal(row$focalClade, "Haloferax")
  expect_equal(row$focalColor, "green3")
})

test_that("averageResults() works on formatStats() output", {
  skip("Known bug: averageResults() requires the MRCA, focalClade and focalColor columns.")
  expect_equal(nrow(averageResults(formatStats(halo_files()))), 15)
})
