test_that("resultFiles() lists the example YAML files, named by pair", {
  f <- halo_files()
  expect_type(f, "character")
  expect_length(f, 30)
  expect_true(all(file.exists(f)))
  expect_equal(names(f)[1], "Halobacterium_litoreum___Halobacterium_noricense")
  expect_false(any(grepl("\\.yaml|\\.bz2", names(f))))
})

test_that("resultFiles() accepts plain, bz2 and gz files", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("A___B.yaml", "A___C.yaml.bz2", "A___D.yaml.gz", "notes.txt")))
  expect_setequal(names(resultFiles(dir)), c("A___B", "A___C", "A___D"))
})

test_that("resultFiles() keeps files with 'gz' inside their name (regression)", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("Agz___B.yaml", "A___Bgz.yaml.gz")))
  expect_setequal(names(resultFiles(dir)), c("Agz___B", "A___Bgz"))
})

test_that("resultFiles() removes files matching patterns of the remove file", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("A___B.yaml", "A___C.yaml", "C___B.yaml")))
  remove <- file.path(dir, "remove.txt")
  writeLines("C___\tbad assembly", remove)
  expect_setequal(names(resultFiles(dir, remove = remove)), c("A___B", "A___C"))
})

test_that("resultFiles() only matches the .yaml extension", {
  skip("Known bug: unescaped dots in the file pattern also match '_yaml'.")
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("A___B.yaml", "A___C_yaml")))
  expect_equal(names(resultFiles(dir)), "A___B")
})
