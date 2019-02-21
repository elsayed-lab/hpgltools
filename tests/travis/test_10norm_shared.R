start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("10norm_shared.R: Are normalizations consistent over time (Shared functions)?\n")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

load("pasilla_df.rda")
test_genes <- c("FBgn0000014","FBgn0000008","FBgn0000017","FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## First make sure the pasilla_expt still has the stuff we expect
expected <- "This is an expt class."
actual <- pasilla_expt[["title"]]
test_that("Pasilla title?", {
    expect_equal(expected, actual)
})

## Ensure that the beginning count table library sizes are identical.
expected <- colSums(counts)
actual <- pasilla_expt[["libsize"]]
names(expected) <- c("untreated1", "untreated2", "untreated3", "untreated4",
                     "treated1", "treated2", "treated3")
test_that("Pasilla libsize?", {
    expect_equal(expected, actual)
})

## Check a few arbitrary counts to make sure they are maintained.
expected <- counts[test_genes, "untreated1"]
testing_counts <- exprs(pasilla_expt)
actual <- as.numeric(testing_counts[test_genes, "untreated1"])
test_that("Pasilla count tables? (untreated1)", {
    expect_equal(expected, actual)
})

## Check that all samples agree for 1 gene.
test_gene <- "FBgn0062565"
expected <- as.numeric(counts[test_gene, ])
actual <- as.numeric(exprs(pasilla_expt)[test_gene, ])
expected <- c(4, 7, 3, 3, 9, 10, 9)
test_that("Pasilla count tables? (gene FBgn0063565)", {
    expect_equal(expected, actual)
})

## Ensure that normalize_expt does not mess up the data when called without arguments (this wasn't true once)
unmolested <- sm(normalize_expt(pasilla_expt))
expected <- as.matrix(exprs(pasilla_expt))  ## I accidently changed this to potentially return a data.frame
actual <- exprs(unmolested)
test_that("Pasilla (un)normalized counts?", {
    expect_equal(expected, actual)
})

## Make sure that the pData information is maintained through normalization
expected <- pData(pasilla_expt)
actual <- pData(unmolested)
test_that("Pasilla (un)normalized pdata?", {
    expect_equal(expected, actual)
})

## Also ensure that the library sizes (which are very important for limma) are not messed up.
expected <- pasilla_expt[["libsize"]]
actual <- unmolested[["libsize"]]
test_that("Pasilla (un)normalized libsize?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 10norm_shared.R in ", elapsed, " seconds."))
tt <- try(clear_session())
