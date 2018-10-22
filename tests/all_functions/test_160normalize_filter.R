start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("160normalize_filter.R:\n")
## 2017-12, exported functions in annotation_gff:

pombe_expt <- make_pombe_expt()

testing <- normalize_expt(pombe_expt, filter=TRUE)
test_counts <- exprs(testing)

expected <- c(23, 37, 155, 19, 91, 184, 49, 105, 151, 22)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does default filtering provide expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})
testing <- normalize_expt(pombe_expt, filter="hpgl")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does the hpgl filter provide expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <- normalize_expt(pombe_expt, filter="pofa")
test_counts <- exprs(testing)

expected <- c(8, 23, 1, 37, 2, 0, 0, 36, 3, 155)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does genefilter's pofa filtering provide expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <- normalize_expt(pombe_expt, filter="kofa")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does genefilter's kofa filtering provide expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <- normalize_expt(pombe_expt, filter="kofa")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
test_that("Does genefilter's cv filtering provide expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 160normalize_filter.R in ", elapsed,  " seconds."))
