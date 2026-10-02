test_that("removeAssemblies() removes matching genomes from matrices", {
  m <- removeAssemblies(Halo_PercentDiff, textConnection("Halobacterium_\ttest"))
  expect_equal(rownames(m), c("Haloferax_mediterranei", "Haloferax_volcanii", "Salarchaeum_japonicum"))
  expect_equal(colnames(m), rownames(m))
})

test_that("removeAssemblies() removes matching pairs from data frames", {
  df <- removeAssemblies(Halo_DF, textConnection("Haloferax_volcanii\tbad\nSalarchaeum\tworse"))
  expect_false(any(grepl("Haloferax_volcanii|Salarchaeum", rownames(df))))
  expect_equal(nrow(df), 12)
})

test_that("removeAssemblies() returns its input when no file is given", {
  expect_identical(removeAssemblies(Halo_DF), Halo_DF)
})

test_that("removeAssemblies() rejects other inputs", {
  expect_error(removeAssemblies(list(a = 1), textConnection("a\tb")), "matrix or a data frame")
})
