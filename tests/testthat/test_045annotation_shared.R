start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("045annotation_shared.R:
  12345\n")
## 2017-12, exported functions in annotation_shared:
## get_gene_size() load_annotations()

## get_gene_size uses load_annotations, and load_annotations just calls one of the others.
## So we will use this as a chance to call two of the more troublesome annotation loaders:
## a orgdb and biomart.

testing <- load_annotations(
    host="useast.ensembl.org",
    type="biomart", species="mmusculus")
expected <- 140000
actual <- nrow(testing[["annotation"]])
## 01
test_that("Do we get some data from load_annotations()?", {
  expect_gt(actual, expected)
})

tmp <- testing[["annotation"]]
testing <- get_genesizes(annotation=tmp, type_column="gene_biotype",
                         gene_type="protein_coding")
actual <- head(testing[["gene_size"]])
expected <- c(1065, 525, 753, 1701, 543, 1038)
## 02
test_that("Do we get consistent gene sizes?", {
  expect_equal(expected, actual)
})

actual <- head(tmp[["start_position"]])
expected <- c(108107280, 77837901, 96271457, 18780447, 142650766, 108343354)
## 03
test_that("Do we get consistent gene starts?", {
  expect_equal(expected, actual)
})

actual <- head(tmp[["end_position"]])
expected <- c(108146146, 77853623, 96276595, 18811987, 142666816, 108414396)
## 04
test_that("Do we get consistent gene ends?", {
  expect_equal(expected, actual)
})

actual <- head(tmp[["transcript_version"]])
expected <- c(4, 13, 8, 13, 11, 5)
## 05
test_that("Do we get consistent transcript versions?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 045annotation_shared.R in ", elapsed,  " seconds."))
