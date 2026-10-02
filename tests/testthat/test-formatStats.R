test_that("formatStats() makes one row per ordered pair of distinct species", {
  df <- formatStats(halo_files())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 30)
  expect_setequal(unique(df$species1), halo_species)
  expect_setequal(unique(df$species2), halo_species)
  expect_false(any(df$species1 == df$species2))
  expect_equal(rownames(df), paste(df$species1, df$species2, sep = "___"))
})

test_that("formatStats() output matches the Halo_DF example data", {
  df <- formatStats(halo_files())
  expect_equal(df, Halo_DF[, names(df)])
})

test_that("formatStats() computes the derived columns", {
  df <- formatStats(halo_files())
  expect_equal(df$percent_difference_global, 100 - df$percent_identity_global)
  expect_equal(df$percent_difference_local,  100 - df$percent_identity_local)
  expect_equal(df$index_avg_strandRand,
               (df$index_strandRand_target + df$index_strandRand_query) / 2)
  expect_equal(df$index_avg_strandDiscord, 1 - df$index_avg_strandRand)
  expect_equal(df$percent_identity_local,
               df$aligned_matches_Total / (df$aligned_matches_Total + df$aligned_mismatches_Total) * 100)
  expect_true(all(df$percent_aligned > 0 & df$percent_aligned <= 100))
})

test_that("formatStats() makes the same label for both orientations of a pair", {
  df <- formatStats(halo_files())
  expect_equal(df["Haloferax_volcanii___Haloferax_mediterranei", "lab"],
               df["Haloferax_mediterranei___Haloferax_volcanii", "lab"])
  expect_equal(df["Haloferax_volcanii___Haloferax_mediterranei", "lab"],
               "Haloferax_volcanii\nHaloferax_mediterranei")
  expect_length(unique(df$lab), 15)
})

test_that("formatStats() reads legacy PercentIdentity fields", {
  stats <- yaml::yaml.load(yaml::read_yaml(halo_files()[1]))
  stats$PercentIdentity       <- 80
  stats$PercentIdentityNoGaps <- 90
  dir <- withr::local_tempdir()
  f <- c(A___B = write_stats_yaml(stats, file.path(dir, "A___B.yaml")))
  df <- formatStats(f)
  expect_equal(df$percent_identity_global, 80)
  expect_equal(df$percent_identity_local,  90)
  expect_equal(df$percent_difference_global, 20)
})

test_that("formatStats() reads legacy PercentSimilarity fields", {
  stats <- yaml::yaml.load(yaml::read_yaml(halo_files()[1]))
  stats$PercentSimilarity       <- 70
  stats$PercentSimilarityNogaps <- 75
  dir <- withr::local_tempdir()
  f <- c(A___B = write_stats_yaml(stats, file.path(dir, "A___B.yaml")))
  df <- formatStats(f)
  expect_equal(df$percent_identity_global, 70)
  expect_equal(df$percent_identity_local,  75)
})

test_that("formatStats() drops self-comparisons", {
  f <- halo_files()[1:2]
  names(f)[2] <- "A___A"
  expect_equal(rownames(formatStats(f)), names(f)[1])
})
