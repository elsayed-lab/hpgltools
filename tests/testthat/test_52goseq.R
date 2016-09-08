library(testthat)
library(hpgltools)

context("Does goseq work?")

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

expected <- c("GO:0004252", "GO:0003824")
actual <- head(rownames(goseq_result$mf_interesting))
test_that("Are the goseq interesting results as expected (mf categories)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0006508")
actual <- head(rownames(goseq_result$bp_interesting))
test_that("Are the goseq interesting results as expected (bp categories)?", {
    expect_equal(expected, actual)
})

expected <- c("GO:0005615", "GO:0016021", "GO:0005576")
actual <- head(rownames(goseq_result$cc_interesting))
test_that("Are the goseq interesting results as expected (cc categories)?", {
    expect_equal(expected, actual)
})

expected <- c(0.2332155, 0.1883853, 0.5333333, 0.2279412, 0.4375000, 0.2129630)
actual <- head(goseq_result$pvalue_plots$mfp_plot_over$data$score)
test_that("Are the goseq results as expected (mf pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.1866913, 0.2135922, 0.1853547, 0.1791444, 0.2075472, 0.2777778)
actual <- head(goseq_result$pvalue_plots$bpp_plot_over$data$score)
test_that("Are the goseq results as expected (bp pvalues)?", {
    expect_equal(expected, actual, tolerance=0.000001)
})

expected <- c(0.1822785, 0.1783088, 0.1606498, 0.1739675, 0.1843972, 0.2471910)
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

## With that in mind, lets test some ontologies:


```{r testing_ontologies_parents}
library(GO.db)
goont("GO:0005576")
golev("GO:0005576")
goont("GO:0048067")
golev("GO:0048067")
goont("GO:0016853")
golev("GO:0016853")
goont("GO:0042438")
golev("GO:0042438")
```
