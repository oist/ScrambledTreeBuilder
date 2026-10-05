## code to prepare `Haloarchaea` data objects goes here

Halo_DF <- system.file("extdata/yaml", package = "ScrambledTreeBuilder") |>
  resultFiles() |>  formatStats()

Halo_PercentDiff <- makeMatrix(Halo_DF, "percent_difference_global", 0)
usethis::use_data(Halo_PercentDiff, overwrite = TRUE)

Halo_ScramblIdx <- makeMatrix(Halo_DF, "index_avg_strandDiscord", 1)

Halo_Tree <- makeTidyTree((Halo_PercentDiff/2 + t(Halo_PercentDiff)/2)) |>
  makeValueTibble(Halo_PercentDiff, colname = "Percent_difference")     |>
  makeValueTibble(Halo_ScramblIdx,  colname = "Scrambling_index")
usethis::use_data(Halo_Tree, overwrite = TRUE)

Halo_bacterium   <- focalClade(Halo_Tree, "Halobacterium_noricense", "Halobacterium_salinarum", "blue",   "Halobacterium")
Halo_ferax       <- focalClade(Halo_Tree, "Haloferax_mediterranei",  "Haloferax_volcanii",      "green3", "Haloferax")
Halo_FocalClades <- FocalCladeList(Halobacterium = Halo_bacterium, Haloferax = Halo_ferax)

Halo_DF <- recordAncestor(Halo_DF, Halo_Tree)
Halo_DF <- recordClades  (Halo_DF, Halo_FocalClades)
usethis::use_data(Halo_DF, overwrite = TRUE)

usethis::use_data(Halo_FocalClades, overwrite = TRUE)

# A made-up species trait, with small integer values like chromosome numbers
# (the real ones are all 1), chosen to illustrate visualizeTree(trait = ...):
# a background value (3), a focal clade with its own value (Haloferax, 2), an
# exception inside the other focal clade (H. litoreum, 4), and a lineage
# (Salarchaeum, 5) whose ancestor parsimony can not resolve.
Halo_Taxons <- data.frame(
  row.names   = c("Halobacterium_litoreum", "Halobacterium_noricense",
                  "Halobacterium_salinarum", "Haloferax_mediterranei",
                  "Haloferax_volcanii", "Salarchaeum_japonicum"),
  Toy_trait = c(4L, 3L, 3L, 2L, 2L, 5L))
usethis::use_data(Halo_Taxons, overwrite = TRUE)
