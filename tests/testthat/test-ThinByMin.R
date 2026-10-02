test_that("thinByMin() keeps labels at least min_dist apart", {
  m <- thinByMin(Halo_PercentDiff, min_dist = 20)
  off_diag <- m[upper.tri(m)]
  expect_true(all(off_diag >= 20))
  expect_equal(rownames(m), colnames(m))
  expect_equal(m, t(m))
  expect_equal(nrow(m), 5)
})

test_that("thinByMin() starts from the most distant pair", {
  m <- thinByMin(Halo_PercentDiff, min_dist = 1000)
  D <- halo_sym()
  expect_equal(nrow(m), 2)
  expect_equal(m[1, 2], max(D))
})

test_that("thinByMin() rejects matrices with NA", {
  m <- Halo_PercentDiff
  m[1, 2] <- NA
  expect_error(thinByMin(m, 20), "NA")
})
