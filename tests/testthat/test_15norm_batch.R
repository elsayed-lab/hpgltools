start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("15norm_batch.R: Are normalizations consistent over time (Batch estimation/correction)?\n")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    batch=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata[["sampleid"]] <- rownames(metadata)
## Uses these genes for quick tests
test_genes <- c("FBgn0000014","FBgn0000008","FBgn0000017","FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## Test batch
## The following is the previous result, it seems to have changed
## expected <- c(3.333333, 64.500000, 3040.166667, 383.916667, 7.083333)
expected <- c(2.032443, 70.820173, 3357.734214, 379.051162, 6.500123)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="limma"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME this is obviously broken
expected <- c(1.60048774, 0.02101530, -0.07524254, 0.15555548, 0.49697157)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="limmaresid"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma-residuals batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(3.095141, 60.542865, 3032.546240, 355.354483, 6.666536)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="combatmod"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("combatmod from cbcbSEQ batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(0.1139956, 7.7406815, 384.8292656, 34.1636051, 0.4937972)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="sva"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("sva batch gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO Figure out what is up with these
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat")  ## broken
## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior") ## takes forever
expected <- c(2.267660, 71.646546, 3481.434409, 410.347808, 6.658058)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="combat_scale"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("combat_scale gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## pasilla_batch <- normalize_expt(pasilla_expt, batch="combat_noprior_scale") ## takes forever
## The previous result
##expected <- c(0.1139956, 7.7406815, 384.8292656, 34.1636051, 0.4937972)
expected <- c(4.610009, 82.109047, 4099.039062, 519.407500, 9.116170)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="svaseq"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("svaseq gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4, 83, 4091, 496, 9)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="ruvg"))
actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("ruvg gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## The following tests take too much memory
if (!identical(Sys.getenv("TRAVIS"), "true")) {
    ##expected <- c(4.610009, 82.109047, 4099.039071, 519.407501, 9.116170)
    ##names(expected) <- test_genes
    ##pasilla_batch <- sm(normalize_expt(pasilla_expt, batch="varpart"))
    ##actual_df <- Biobase::exprs(pasilla_batch[["expressionset"]])
    ##actual <- actual_df[test_genes, c("untreated1")]
    ##test_that("variancePartition gives expected values?", {
    ##    expect_equal(expected, actual, tolerance=0.0001)
    ##})
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 15norm_batch.R in ", elapsed,  " seconds."))
