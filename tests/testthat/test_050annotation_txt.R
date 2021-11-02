start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("050annotation_txt.R:
  123\n")
## 2017-12, exported functions in annotation_txt:
## load_trinotate_annotations() load_trinotate_go()

tmp <- system.file("share/sb/trinotate_head.csv.xz", package = "hpgltools")
testing <- load_trinotate_annotations(trinotate = tmp)
## Moved rownames to a column for tibble, ergo 1 more column
## I recently changed this to drop columns which have no information, losing 3.
expected <- c(49999, 31)
actual <- dim(testing)
## 0102
test_that("Do we get expected trinotate data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## FIXME: The creation of this go table is not quite right I think.
testing <- load_trinotate_go(trinotate = tmp)
## expecteed <- c(136152, 5)
expected <- c(112023, 5)
actual <- dim(testing[["go_data"]])
## 03
test_that("Do we get expected trinotate GO data?", {
  expect_equal(expected[1], actual[1])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 050annotation_txt.R in ", elapsed,  " seconds.")
