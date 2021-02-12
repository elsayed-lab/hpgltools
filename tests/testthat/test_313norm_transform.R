start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("13norm_transform.R: Are normalizations consistent over time (Transformations)?\n")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

load("pasilla_df.rda")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]

## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.

## Test transformations
expected <- c(2.584963, 6.539159, 12.187661, 9.189825, 3.459432)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform = "log2"))
actual_df <- exprs(pasilla_trans)
actual <- actual_df[test_genes, c("untreated1")]
test_that("log2 transformation gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(0.7781513, 1.9684829, 3.6688516, 2.7664128, 1.0413927)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform = "log10"))
actual_df <- exprs(pasilla_trans)
actual <- actual_df[test_genes, c("untreated1")]
test_that("log10 transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(1.791759, 4.532599, 8.447843, 6.369901, 2.397895)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform = "log"))
actual_df <- exprs(pasilla_trans)
actual <- actual_df[test_genes, c("untreated1")]
test_that("loge transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 13norm_transform.R in ", elapsed,  " seconds."))
