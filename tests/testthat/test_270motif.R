start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("270motif.R")

## Also going to fail because eupathdb is not working yet.

if (FALSE) {
  tritryp_metadata <- EuPathDB::download_eupath_metadata(webservice = "tritrypdb")
  lm_entry <- EuPathDB::get_eupath_entry(species = "major", metadata = tritryp_metadata)
  lm_org <- EuPathDB::make_eupath_bsgenome(entry = lm_entry)
  pkgnames <- EuPathDB::get_eupath_pkgnames(entry = lm_entry)
  lm_bsg <- pkgnames[["bsgenome"]]

  motif_input <- system.file("share/motif/pro_high.fasta", package = "hpgldata")

  lm_gadem <- sm(simple_gadem(motif_input, genome = lm_bsg, p = 1.0))

  actual <- lm_gadem[["pvals"]][[1]]
  expected <- 5e-5
  test_that("Do we get an expected p-value?", {
    expect_lt(actual, expected)
  })

}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 270motif.R in ", elapsed,  " seconds.")
