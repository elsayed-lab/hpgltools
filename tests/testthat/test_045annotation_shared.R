start <- as.POSIXlt(Sys.time())
context("045annotation_shared.R")
## 2017-12, exported functions in annotation_shared:
## get_gene_size() load_annotations()

## get_gene_size uses load_annotations, and load_annotations just calls one of the others.

pa_gff <- system.file("share/paeruginosa_pa14.gff", package = "hpgldata")
pa_fasta <- system.file("share/paeruginosa_pa14.fasta", package = "hpgldata")

testing <- load_annotations(type = "gff", gff = pa_gff)
expected <- 11946
actual <- nrow(testing)
## 01
test_that("Do we get some data from load_annotations()?", {
  expect_equal(actual, expected)
})

testing <- get_genesizes(annotation = testing)
actual <- head(testing[["gene_size"]])
expected <- c(1545, 1104, 1110, 2421, 774, 537)
## 02
test_that("Do we get consistent gene sizes?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 045annotation_shared.R in ", elapsed,  " seconds.")
