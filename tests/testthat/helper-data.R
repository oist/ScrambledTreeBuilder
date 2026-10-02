# Shared fixtures.  The lazy-loaded Halo_* objects are rebuilt by
# data-raw/Haloarchaea.R from the YAML files in inst/extdata/yaml.

halo_files <- function() {
  system.file("extdata/yaml", package = "ScrambledTreeBuilder") |> resultFiles()
}

halo_species <- c("Halobacterium_litoreum", "Halobacterium_noricense",
                  "Halobacterium_salinarum", "Haloferax_mediterranei",
                  "Haloferax_volcanii", "Salarchaeum_japonicum")

# Symmetrised percent difference matrix, as used to build Halo_Tree.
halo_sym <- function() (Halo_PercentDiff + t(Halo_PercentDiff)) / 2

# Write a list of statistics as a YAML file in the format produced by the
# pipeline: a YAML document containing a single string that is itself YAML.
write_stats_yaml <- function(stats, path) {
  yaml::write_yaml(yaml::as.yaml(stats), path)
  path
}

# Pretend chromosome numbers for the Haloarchaea (the real ones are all 1).
# Haloferax (2) is the outgroup of Salarchaeum (5) + Halobacterium (3), so
# parsimony cannot resolve node 9 (3 or 5) nor the root (2, 3 or 5).
halo_taxons <- function() {
  data.frame(row.names = halo_species,
             Binomial    = sub("_", " ", halo_species),
             ChromNumber = c(3, 3, 3, 2, 2, 5))
}

# A four-tip tree ((A,B),(C,D)), built from a distance matrix.
abcd_tree <- function() {
  m <- matrix(c( 0,  2, 10, 10,
                 2,  0, 10, 10,
                10, 10,  0,  4,
                10, 10,  4,  0), 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  makeTidyTree(m)
}

# Node ID of the most recent common ancestor of two tips.
mrca_of <- function(tree, a, b) {
  treeio::MRCA(tree, which(tree$label %in% a), which(tree$label %in% b))$node
}
