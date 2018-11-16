start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("005annotation_biomart.R
  1234\n")
## 2017-12, exported functions in annotation_biomart:
##   load_biomart_annotations(), load_biomart_go(), load_biomart_orthologs()

## load_biomart_annotations()
testing <- load_biomart_annotations(species="hsapiens")
gene_ids <- head(rownames(testing[["annotation"]]))
expected <- c("ENST00000000233", "ENST00000000412", "ENST00000000442",
              "ENST00000001008", "ENST00000001146", "ENST00000002125")
## 01
test_that("Do we get expected gene IDs?", {
  expect_equal(expected, gene_ids)
})

data <- testing[["annotation"]]
expected <- c(197995, 12)
actual <- dim(data)
## 02
test_that("Do we receive expected output from load_biomart_annotations()?", {
  expect_equal(expected, actual)
})

## load_biomart_go()
testing <- load_biomart_go(species="hsapiens")
data <- testing[["go"]]
expected <- c(318558, 2)
actual <- dim(data)
## 03
test_that("Do we receive expected output from load_biomart_go()?", {
  expect_equal(expected, actual)
})

## load_biomart_orthologs()
## Oh yeah, I moved the default biomart to hg38/90 or 89.  Thus these are unlikely
## to be correct.  I should just query the number of orthologs found.
testing <- load_biomart_orthologs(gene_ids=gene_ids, first_species="hsapiens",
                                  second_species="mmusculus")
data <- testing[["all_linked_genes"]]
actual <- dim(data)
expected <- c(26635, 3)
## 04
test_that("Do we get expected orthologs from load_biomart_orthologs()?", {
  expect_equal(expected[1], actual[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 005annotation_biomart.R in ", elapsed,  " seconds."))
