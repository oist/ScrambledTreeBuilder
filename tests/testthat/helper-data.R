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
