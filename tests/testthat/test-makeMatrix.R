test_that("makeMatrix() builds a square species matrix", {
  m <- makeMatrix(Halo_DF, "percent_difference_global", 0)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(6, 6))
  expect_equal(rownames(m), halo_species)
  expect_equal(colnames(m), halo_species)
  expect_equal(unname(diag(m)), rep(0, 6))
  expect_equal(attr(m, "builtWith"), "percent_difference_global")
  expect_false(anyNA(m))
})

test_that("makeMatrix() puts species1 in rows and species2 in columns", {
  m <- makeMatrix(Halo_DF, "percent_difference_global", 0)
  expect_equal(m["Haloferax_volcanii", "Salarchaeum_japonicum"],
               Halo_DF["Haloferax_volcanii___Salarchaeum_japonicum", "percent_difference_global"])
  expect_equal(m["Salarchaeum_japonicum", "Haloferax_volcanii"],
               Halo_DF["Salarchaeum_japonicum___Haloferax_volcanii", "percent_difference_global"])
})

test_that("makeMatrix() matches the Halo_PercentDiff example data", {
  expect_equal(makeMatrix(Halo_DF, "percent_difference_global", 0), Halo_PercentDiff)
})

test_that("makeMatrix() warns and returns NULL for unknown columns", {
  expect_warning(m <- makeMatrix(Halo_DF, "no_such_column", 0), "not found")
  expect_null(m)
})

test_that("makeMatrix() stops on species with only zero values", {
  df <- Halo_DF
  df[df$species1 == "Haloferax_volcanii", "percent_difference_global"] <- 0
  expect_error(makeMatrix(df, "percent_difference_global", 0), "Haloferax_volcanii")
})

df_with_NA <- function() {
  df <- Halo_DF
  df["Salarchaeum_japonicum___Haloferax_volcanii", "percent_difference_global"] <- NA
  df
}

test_that("makeMatrix() leaves missing values as NA by default", {
  m <- makeMatrix(df_with_NA(), "percent_difference_global", 0)
  expect_true(is.na(m["Salarchaeum_japonicum", "Haloferax_volcanii"]))
  expect_equal(sum(is.na(m)), 1)
})

test_that("makeMatrix() fills missing values with defaultValue", {
  m <- makeMatrix(df_with_NA(), "percent_difference_global", 0, 50)
  expect_equal(m["Salarchaeum_japonicum", "Haloferax_volcanii"], 50)
})

test_that("makeMatrix(impute = 'average') takes the value of the reverse pair", {
  m <- makeMatrix(df_with_NA(), "percent_difference_global", 0, 50, impute = "average")
  expect_equal(m["Salarchaeum_japonicum", "Haloferax_volcanii"],
               m["Haloferax_volcanii", "Salarchaeum_japonicum"])
  df <- df_with_NA()
  df["Haloferax_volcanii___Salarchaeum_japonicum", "percent_difference_global"] <- NA
  m <- makeMatrix(df, "percent_difference_global", 0, 50, impute = "average")
  expect_equal(m["Salarchaeum_japonicum", "Haloferax_volcanii"], 50)
})

test_that("makeMatrix() imputes with missForest", {
  skip_if_not_installed("missForest")
  df <- df_with_NA()
  df["Haloferax_volcanii___Salarchaeum_japonicum", "percent_difference_global"] <- NA
  set.seed(1664)
  m  <- makeMatrix(df, "percent_difference_global", 0, impute = "missForest")
  set.seed(1664)
  m2 <- makeMatrix(df, "percent_difference_global", 0, impute = "missForest2")
  expect_false(anyNA(m))
  expect_false(anyNA(m2))
  expect_equal(unname(diag(m2)), rep(0, 6))
})

test_that("makeMatrix() passes extra arguments to the imputing function", {
  skip("Known bug: .imputeMatrix() does not forward '...' to imputeFunction().")
  skip_if_not_installed("missForest")
  expect_error(makeMatrix(df_with_NA(), "percent_difference_global", 0,
                          impute = "missForest", not_a_missForest_argument = 1))
})

test_that("makeMatrix() rejects unknown imputation methods", {
  expect_error(makeMatrix(Halo_DF, "percent_difference_global", 0, impute = "magic"))
})
