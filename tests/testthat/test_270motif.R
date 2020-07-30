start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("270motif.R:
")

tritryp_metadata <- EuPathDB::download_eupath_metadata(webservice="tritrypdb")
lm_entry <- EuPathDB::get_eupath_entry(species="major", metadata=tritryp_metadata)
lm_org <- EuPathDB::make_eupath_bsgenome(entry=lm_entry)
pkgnames <- EuPathDB::get_eupath_pkgnames(entry=lm_entry)
lm_bsg <- pkgnames[["bsgenome"]]

motif_input <- system.file("motif/pro_high.fasta", package="hpgltools")

lm_gadem <- simple_gadem(motif_input, genome=lm_bsg, p=1.0)

actual <- lm_gadem[["pvals"]][["CGTGCGTGTG"]]
expected <- 4.686288e-05
test_that("Do we get an expected p-value?", {
  expect_equal(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 270motif.R in ", elapsed,  " seconds."))
