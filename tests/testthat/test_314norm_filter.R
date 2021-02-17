start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("14norm_filter.R: Are normalizations consistent over time (Filtering)?\n")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.
## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.

load("pasilla_df.rda")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]

## Test filter
expected <- c(7531, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter = "cbcb"))
actual <- dim(exprs(pasilla_filter))
test_that("cbcb filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter = "pofa"))
actual <- dim(exprs(pasilla_filter))
test_that("pofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter = "kofa"))
actual <- dim(exprs(pasilla_filter))
test_that("kofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter = "cv"))
actual <- dim(exprs(pasilla_filter))
test_that("cv filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter = "simple"))
actual <- dim(exprs(pasilla_filter))
test_that("simple filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 14norm_filter.R in ", elapsed,  " seconds.")
