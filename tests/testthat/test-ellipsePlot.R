test_that("ellipsePlot() plots pairs grouped by MRCA and colored by clade", {
  a <- averageResults(Halo_DF)
  p <- ellipsePlot(a)
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2)
  d <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(d), 15)
  expect_equal(d$x, a$percent_difference_local)
  expect_true(all(d$colour[a$focalClade %in% "Haloferax"] == "green3"))
})

test_that("ellipsePlot() uses the requested statistics as axis labels", {
  p <- ellipsePlot(averageResults(Halo_DF), x = "percent_difference_global", y = "percent_aligned")
  expect_equal(p$labels$x, "percent_difference_global")
  expect_equal(p$labels$y, "percent_aligned")
})
