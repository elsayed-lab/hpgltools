start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)

context("52gsea_goseq.R: Does goseq work?\n")
## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.

##if (!identical(Sys.getenv("TRAVIS"), "true")) {
load("gsea_siggenes.rda")

goseq_result <- sm(simple_goseq(de_genes=fcp_sig_genes, length_db=dmel_lengths, go_db=dmel_ontologies))

expected <- 53
actual <- nrow(goseq_result$mf_interesting)
test_that("Do we get the expected number of interesting mf categories?", {
    expect_equal(expected, actual, tolerance=2)
})
expected <- 71
actual <- nrow(goseq_result$bp_interesting)
test_that("Do we get the expected number of interesting bp categories?", {
    expect_equal(expected, actual)
})
expected <- 13
actual <- nrow(goseq_result$cc_interesting)
test_that("Do we get the expected number of interesting cc categories?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0000146", "GO:0003824", "GO:0003974",
              "GO:0003978", "GO:0003993", "GO:0004114")
actual <- head(sort(rownames(goseq_result$mf_interesting)))
test_that("Are the goseq interesting results as expected (mf categories)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0000422", "GO:0000578", "GO:0001508",
              "GO:0001676", "GO:0002121", "GO:0002385")
actual <- head(sort(rownames(goseq_result$bp_interesting)))
test_that("Are the goseq interesting results as expected (bp categories)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0005576", "GO:0005578", "GO:0005639",
              "GO:0005811", "GO:0005859", "GO:0005967")
actual <- head(sort(rownames(goseq_result$cc_interesting)))
test_that("Are the goseq interesting results as expected (cc categories)?", {
    expect_equal(expected, actual)
})

expected <- c(0.02857143, 0.02884615, 0.02941176, 0.03333333, 0.04651163, 0.07142857)
actual <- head(goseq_result$pvalue_plots$mfp_plot_over$data$score)
test_that("Are the goseq results as expected (mf pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.07692308, 0.08333333, 0.09090909, 0.09090909, 0.10000000, 0.10000000)
actual <- head(goseq_result$pvalue_plots$bpp_plot_over$data$score)
test_that("Are the goseq results as expected (bp pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.01263538, 0.01895735, 0.03636364, 0.03846154, 0.09090909, 0.11111111)
actual <- head(goseq_result$pvalue_plots$ccp_plot_over$data$score)
test_that("Are the goseq results as expected (cc pvalues)?", {
    expect_equal(expected, actual, tolerance=0.01)
})

##}

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
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 52gsea_goseq.R in ", elapsed,  " seconds."))
