start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("12norm_normalize.R: Are normalizations consistent over time? (Normalizations)?
  1234567\n")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

load("pasilla_df.rda")
## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

expected <- as.numeric(c(5.857143, 91.500000, 4400.000000, 543.785714, 10.714286))
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="quant"))
actual_df <- exprs(pasilla_norm)
actual <- as.numeric(actual_df[test_genes, c("untreated1")])
test_that("quant normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=20)
})

## Similar test for size-factor normalization
expected <- c(4.392658, 80.824908, 4097.471407, 512.183926, 8.785316)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="sf"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("size-factor normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Check another size-factor normalization
expected <- c(4.392658, 80.824908, 4097.471407, 512.183926, 8.785316)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="sf2"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("size-factor2 normalization gives expected values?", {
    expect_equal(expected, actual)
})
## Oh I never noticed before that this is a log, too
## Updated this to use DESeq2 estimateDispersions.
##expected <- c(5.488150, 7.082043, 12.021996, 9.160395, 5.707992)
expected <- c(6.390378, 7.554126, 12.045030, 9.310557, 6.544189)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="vsd"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("vsd normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.927997, 91.830657, 4765.366532, 613.466245, 9.342734)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="tmm"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("tmm normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.902692, 90.336774, 4803.090308, 608.726226, 9.488822)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="upperquartile"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("upperquartile normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.927854, 91.079703, 4840.296148, 615.582521, 9.205998)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="rle"))
actual_df <- exprs(pasilla_norm)
actual <- actual_df[test_genes, c("untreated1")]
test_that("RLE normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 12norm_normalize.R in ", elapsed, " seconds."))
tt <- try(clear_session())
