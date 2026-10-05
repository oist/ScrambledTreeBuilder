bar_layer <- function(p) which(sapply(p$layers, \(l) class(l$geom)[1] == "GeomRect"))

test_that("cladeBars() draws one bar per clade next to its tips", {
  p <- visualizeTree(Halo_Tree) + cladeBars(Halo_FocalClades)
  expect_s3_class(p, "ggtree")
  expect_length(bar_layer(p), 1)
  d <- ggplot2::layer_data(p, bar_layer(p))
  expect_equal(nrow(d), 2)
  tips <- p$data[p$data$isTip, ]
  for (i in 1:2) {
    y <- tips$y[tips$label %in% Halo_FocalClades[[i]]@genomeIDs]
    expect_equal(d$ymin[i], min(y) - 0.4)
    expect_equal(d$ymax[i], max(y) + 0.4)
  }
  expect_true(all(d$xmin > max(p$data$x)))
  expect_true(all(d$xmax > d$xmin))
})

test_that("cladeBars() uses the clade colors and names", {
  p <- visualizeTree(Halo_Tree) + cladeBars(Halo_FocalClades)
  d <- ggplot2::layer_data(p, bar_layer(p))
  expect_equal(d$fill, c("blue", "green3"))
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
  expect_equal(scale$get_labels(), c("Halobacterium", "Haloferax"))
})

test_that("cladeBars() accepts a single focal clade", {
  p <- visualizeTree(Halo_Tree) + cladeBars(Halo_FocalClades$Haloferax)
  expect_equal(nrow(ggplot2::layer_data(p, bar_layer(p))), 1)
})

test_that("cladeBars() skips clades absent from the plotted tree", {
  p <- visualizeTree(subTree(Halo_Tree, 9)) + cladeBars(Halo_FocalClades)
  d <- ggplot2::layer_data(p, bar_layer(p))
  expect_equal(nrow(d), 1)
  expect_equal(d$fill, "blue")
})

test_that("cladeBars() are placed after the karyotype number column", {
  p <- visualizeTree(halo_chr_tree(), value = NULL, trait = "ChromNumber") + cladeBars(Halo_FocalClades)
  numbers <- ggplot2::layer_data(p, which(sapply(p$layers, \(l) class(l$geom)[1] == "GeomText")))
  bars <- ggplot2::layer_data(p, bar_layer(p))
  expect_true(all(bars$xmin > max(numbers$x)))
  # The fill scale of the bars does not disturb the colour scale of the branches.
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(scale$get_labels(), c("2", "3", "5", "ambiguous"))
})

test_that("cladeBars() position and width can be set", {
  p <- visualizeTree(Halo_Tree) + cladeBars(Halo_FocalClades, offset = 0.5, width = 0.1)
  d <- ggplot2::layer_data(p, bar_layer(p))
  xmax <- max(p$data$x)
  expect_equal(d$xmin, rep(xmax * 1.5, 2))
  expect_equal(d$xmax - d$xmin, rep(xmax * 0.1, 2))
})

test_that("cladeBars() puts nested clades in separate columns", {
  clades <- FocalCladeList(
    Halobacteria = focalClade(Halo_Tree, "Salarchaeum_japonicum", "Halobacterium_salinarum", "grey", "Halobacteria"),
    Halobacterium = Halo_FocalClades$Halobacterium,
    Haloferax     = Halo_FocalClades$Haloferax)
  p <- visualizeTree(Halo_Tree) + cladeBars(clades)
  d <- ggplot2::layer_data(p, bar_layer(p))
  # The largest clade is closest to the tree; the nested one is next to it.
  expect_equal(d$xmin[1], d$xmin[3])
  expect_gt(d$xmin[2], d$xmax[1])
  # No two bars in the same column overlap.
  for (i in 1:2) for (j in (i + 1):3)
    if (d$xmin[i] == d$xmin[j])
      expect_true(d$ymax[i] < d$ymin[j] || d$ymax[j] < d$ymin[i])
})
