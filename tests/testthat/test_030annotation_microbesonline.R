start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("030annotation_microbesonline.R
  1234\n")
## 2017-12, exported functions in annotation_microbesonline:
## get_microbesonline_ids(), get_microbesonline_name(), load_microbesonline_annotations()
## load_microbesonline_go(), load_microbesonline_kegg().

## ## load_microbesonline_annotations()
spy_annot <- load_microbesonline_annotations(id = 160490)
actual <- dim(spy_annot)
expected <- c(1871, 18)
## 0102
test_that("Do we get the expected amount of pyogenes data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## ## load_microbesonline_go()
spy_go <- load_microbesonline_go(id = 160490)
actual <- dim(spy_go)
expected <- c(4594, 2)
## 0304
test_that("Do we get the expected amount of ecoli GO data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## ## load_microbesonline_kegg()
end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 030nnotation_microbesonline.R in ", elapsed,  " seconds.")
