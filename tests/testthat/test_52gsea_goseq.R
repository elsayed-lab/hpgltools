library(testthat)
library(hpgltools)

context("52gsea_goseq.R: Does goseq work?\n")

## I want to do some much longer tests using goseq/clusterprofiler/topgo/gostats/gprofiler
## These will likely not work with Travis as they take forever.
## I am not however, certain how to use skip_on_travis(), so I printed it and am copying the
## useful bits here.


##if (!identical(Sys.getenv("TRAVIS"), "true")) {
limma <- new.env()
load("de_limma.rda", envir=limma)
table <- limma$hpgl_table
sig_genes <- sm(get_sig_genes(table, column="untreated")$up_genes)

dmel_annotations <- sm(get_biomart_annotations(species="dmelanogaster"))
dmel_ontologies <- sm(get_biomart_ontologies(species="dmelanogaster"))
dmel_lengths <- dmel_annotations[, c("geneID", "length")]
colnames(dmel_lengths) <- c("ID","width")
rownames(dmel_lengths) <- make.names(dmel_lengths[["ID"]], unique=TRUE)
## Drop all duplicate gene IDs
dmel_lengths <- dmel_lengths[ !grepl("\\.", rownames(dmel_lengths)), ]
goseq_result <- sm(simple_goseq(de_genes=sig_genes, length_db=dmel_lengths, go_db=dmel_ontologies))

expected <- c("GO:0003824")
actual <- head(rownames(goseq_result$mf_interesting))
test_that("Are the goseq interesting results as expected (mf categories)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0008152", "GO:0048096")
actual <- head(rownames(goseq_result$bp_interesting))
test_that("Are the goseq interesting results as expected (bp categories)?", {
    expect_equal(expected, actual)
})

##expected <- c("GO:0005615", "GO:0016021", "GO:0005576")
##actual <- head(rownames(goseq_result$cc_interesting))
##test_that("Are the goseq interesting results as expected (cc categories)?", {
##    expect_equal(expected, actual)
##})

expected <- c(0.1430595, 0.3750000, 0.1732283, 0.1307420, 0.2280702, 0.2000000)
actual <- head(goseq_result$pvalue_plots$mfp_plot_over$data$score)
test_that("Are the goseq results as expected (mf pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.1684492, 0.1533181, 0.2777778, 0.2692308, 0.4615385, 0.4117647)
actual <- head(goseq_result$pvalue_plots$bpp_plot_over$data$score)
test_that("Are the goseq results as expected (bp pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.2340426, 0.3125000, 0.1421801, 0.1797753, 0.3636364, 0.1666667)
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

message("\nFinished 52gsea_goseq.R")
