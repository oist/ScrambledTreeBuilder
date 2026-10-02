test_that("childSpecies() lists the tips under each branch of a node", {
  expect_equal(childSpecies(Halo_Tree, 8),
               list(left = "Haloferax_mediterranei", right = "Haloferax_volcanii"))
  root <- childSpecies(Halo_Tree, 7)
  expect_setequal(root$left,  c("Haloferax_mediterranei", "Haloferax_volcanii"))
  expect_setequal(root$right, setdiff(halo_species, root$left))
})

test_that("childSpecies() fails on tips", {
  expect_error(childSpecies(Halo_Tree, 1))
})
