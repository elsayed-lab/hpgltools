start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("52gsea_goseq.R: Does goseq work?\n")
load("gsea_siggenes.rda")

goseq_result <- simple_goseq(fcp_sig_genes,
                             length_db=dmel_lengths,
                             go_db=dmel_ontologies)

expected <- 53
actual <- nrow(goseq_result[["mf_interesting"]])
test_that("Do we get the expected number of interesting mf categories?", {
    expect_equal(expected, actual, tolerance=2)
})

expected <- 71
actual <- nrow(goseq_result[["bp_interesting"]])
test_that("Do we get the expected number of interesting bp categories?", {
    expect_equal(expected, actual, tolerance=2)
})

expected <- 8
actual <- nrow(goseq_result[["cc_interesting"]])
test_that("Do we get the expected number of interesting cc categories?", {
    expect_equal(expected, actual)
})

expected <- 45
actual <- nrow(goseq_result[["mf_interesting"]])
test_that("Are the goseq interesting results as expected (mf categories)?", {
    expect_gt(actual, expected)
})

expected <- 60
actual <- nrow(goseq_result[["bp_interesting"]])
test_that("Are the goseq interesting results as expected (bp categories)?", {
    expect_gt(actual, expected)
})

expected <- 5
actual <- nrow(goseq_result[["cc_interesting"]])
test_that("Are the goseq interesting results as expected (cc categories)?", {
    expect_gt(actual, expected)
})

## New goseq versions get slightly lower p-values
## ^^ That was untrue, it was actually due to slightly different numbers of genes.
expected <- 0.14
summ <- summary(goseq_result[["pvalue_plots"]][["mfp_plot_over"]][["data"]][["score"]])
actual <- summ[["Mean"]]
test_that("Are the goseq results as expected (mf pvalues)?", {
    expect_gt(expected, actual)
})

actual <- summary(goseq_result[["pvalue_plots"]][["bpp_plot_over"]][["data"]][["score"]])
actual <- actual[["Mean"]]
test_that("Are the goseq results as expected (bp pvalues)?", {
    expect_gt(expected, actual)
})

## Some testing of an interesting point by keith:
## I've been trying to figure out why goseq is giving me a different gene/term mapping that the one I create myself.
## I think I finally tracked it down to this:
##    #Because GO is a directed graph, we need to get not just the genes associated with each ID,
##    #but also those associated with its children.  GO2ALLEGS does this.
## in: https://github.com/Bioconductor-mirror/goseq/blob/master/R/getgo.R
##
## I know you aren't really using goseq much these days, but I just wanted to run this by you.
## Does this seem right to you?
## The result is that goseq associates a bunch of genes with each go term which aren't actually directly associated with that term.
## So when I think there are 25 genes that have the annotation "immune response", goseq will think there are 75.
## This was not my understanding of how enrichment analysis works. Have I just been mistaken this whole time?
## Keith

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 52gsea_goseq.R in ", elapsed,  " seconds."))
