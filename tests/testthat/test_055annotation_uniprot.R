start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("055annotation_uniprot.R:
  1234\n")
## 2017-12, exported functions in annotation_gff:
## load_uniprot_annotations() download_uniprot_proteome()

## load_uniprot_annotations()
testing <- load_uniprot_annotations(species = "H37Rv")
expected <- c(3993, 76)
actual <- dim(testing)
## 0304
test_that("Did we get the correct uniprot annotations?", {
  expect_gt(actual[1], expected[1])
  expect_gt(actual[2], expected[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 055annotation_uniprot.R in ", elapsed,  " seconds.")
