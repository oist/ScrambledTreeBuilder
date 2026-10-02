test_that("remove_outliers() hides both orientations of outlier pairs", {
  r <- remove_outliers(Halo_DF, "percent_identity_local", "percent_aligned", 72, 40)
  expect_equal(rownames(r)[is.na(r$percent_identity_local)],
               c("Haloferax_volcanii___Salarchaeum_japonicum",
                 "Salarchaeum_japonicum___Haloferax_volcanii"))
  expect_equal(r[, names(r) != "percent_identity_local"],
               Halo_DF[, names(Halo_DF) != "percent_identity_local"])
})

test_that("remove_outliers() returns the input when there are no outliers", {
  expect_equal(remove_outliers(Halo_DF, "percent_identity_local", "percent_aligned", 100, 0), Halo_DF)
})

test_that("remove_outliers() computes default thresholds from quantiles", {
  r <- remove_outliers(Halo_DF, "percent_identity_local", "percent_aligned")
  expect_equal(nrow(r), nrow(Halo_DF))
})
