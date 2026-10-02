# Synthetic MRCA summaries along a noisy linear trend, in four clades of ten
# points.  Clade "A" (low x) is shifted above the trend.
enr_points <- function() {
  set.seed(3)
  tb <- tibble::tibble(x = 1:40, y = 0.02 * (1:40) + rnorm(40, sd = 0.02),
                       n = 2L, clade = rep(c("A", "B", "C", "D"), each = 10))
  tb$y[3:6] <- tb$y[3:6] + 0.15
  tb
}

# Same, but with clade "B" below the trend and clade "D" (high x) above it.
enr_points_2 <- function() {
  tb <- enr_points()
  tb$y[3:6]   <- tb$y[3:6]   - 0.15
  tb$y[13:16] <- tb$y[13:16] - 0.15
  tb$y[33:36] <- tb$y[33:36] + 0.15
  tb
}

test_that("computeENR() adds fitted values, residuals and z-scores", {
  tb <- computeENR(enr_points())
  expect_named(tb, c(names(enr_points()), "y_hat", "resid", "ENRz"))
  expect_equal(tb$resid, tb$y - tb$y_hat)
  expect_true(all(tb$ENRz[3:6] > 3))
  expect_true(all(abs(tb$ENRz[-(3:6)]) < 1.5))
})

test_that("computeENR() standardisation modes differ only in scale", {
  for (mode in c("robust_global", "global_sd", "local_loess")) {
    tb <- computeENR(enr_points_2(), std_mode = mode)
    expect_false(anyNA(tb$ENRz), label = mode)
    expect_equal(sign(tb$ENRz), sign(tb$resid), label = mode)
    expect_lt(tb$ENRz[14], -1)
    expect_gt(tb$ENRz[34],  1)
  }
  expect_error(computeENR(enr_points(), std_mode = "nope"))
})

test_that("computeENR() requires x and y", {
  expect_error(computeENR(data.frame(a = 1)))
})

test_that("cladeENRtable() labels clades above the trend", {
  tbl <- cladeENRtable(computeENR(enr_points()))
  expect_equal(tbl$clade, c("A", "D", "B", "C"))
  expect_equal(tbl$ENR_label, c("ABOVE-TREND", "TYPICAL", "TYPICAL", "TYPICAL"))
  expect_equal(tbl$n_nodes, rep(20L, 4))
  expect_equal(tbl$x_mean[tbl$clade == "A"], 5.5)
})

test_that("cladeENRtable() labels clades below the trend", {
  tbl <- cladeENRtable(computeENR(enr_points_2()))
  expect_equal(tbl$ENR_label[tbl$clade == "B"], "BELOW-TREND")
})

test_that("cladeENRtable() requires above-trend clades to have low divergence", {
  tbl <- cladeENRtable(computeENR(enr_points_2()))
  expect_equal(tbl$ENR_label[tbl$clade == "D"], "TYPICAL")
  tbl <- cladeENRtable(computeENR(enr_points_2()), require_low_x = FALSE)
  expect_equal(tbl$ENR_label[tbl$clade == "D"], "ABOVE-TREND")
})

test_that("cladeENRtable() ignores the 'Other' clade", {
  pts <- enr_points()
  pts$clade[pts$clade == "C"] <- "Other"
  expect_false("Other" %in% cladeENRtable(computeENR(pts))$clade)
})

test_that("cladeENRtable() requires ENR columns", {
  expect_error(cladeENRtable(enr_points()), "must contain")
})

test_that("computeENR() works on MRCAs() output", {
  skip("Known bug: MRCAs() now includes pair rows with n = NA, which break the LOESS weights.")
  tb <- ScrambledTreeBuilder:::MRCAs(Halo_DF, Halo_FocalClades)
  expect_no_error(computeENR(tb))
})
