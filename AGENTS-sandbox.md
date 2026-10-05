# Notes for sandboxed sessions

Read this if the R packages are missing, for instance in a bubblewrap
sandbox where the home directory and the session history may be lost.

## Installing the R dependencies

Install into the user library (`Sys.getenv("R_LIBS_USER")`), not
system-wide.  On Debian with R 4.6 (2026):

- The system BiocManager is pinned to an old Bioconductor and refuses to
  switch without upgrading the system packages.  Bypass it:
  `install.packages(pkgs, lib = Sys.getenv("R_LIBS_USER"), repos = c(CRAN =
  "https://cloud.r-project.org", BioCsoft =
  "https://bioconductor.org/packages/3.23/bioc"))`.
- The system `cpp11` is too old for R 4.6: install `cpp11` first, then
  `roxygen2` and `devtools`.
- Packages that were missing: ggnewscale, plotly, tidytree, missForest,
  patchwork, ggtree and treeio (Bioconductor), roxygen2, devtools.

## Checking

- `devtools::build_rmd()` hung for over 20 minutes.  To check the
  vignette, install to a temporary library and render it there:
  `R CMD INSTALL -l $TMP/lib .`, then
  `R_LIBS=$TMP/lib Rscript -e 'rmarkdown::render("vignettes/ScrambledTreeBuilder.Rmd")'`
  from a copy outside the repository.
- `devtools::run_examples()` leaves an `Rplots.pdf` in the package root:
  delete it.
- Write scratch files outside the repository.

## Real data

`Ovalentaria.RData` and `Eupercaria.RData` (one per branch of the
Neoteleostei) contain: `Tibble` (the `stbTree` built from the `P`
distance, with values from all matrices), `taxons` (taxon table with
`ChromNumber`, `Family`, …), `clades` (a `FocalCladeList`),
`matrices_symmetric` (an `S4Vectors::SimpleList` of species × species
matrices; load S4Vectors) and the pairwise results (`resultsSym`, or
`results` with both directions of each pair in Eupercaria).  Example:
`subsampleSpecies(Tibble, "P", replicates = 5, breaks = seq(0, 0.24,
0.02), clades = clades)`.
