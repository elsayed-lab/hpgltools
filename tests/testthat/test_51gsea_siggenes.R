start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("51gsea_siggenes.R: Do we get consistent sets of 'significant' genes for GSEA analyses?\n")

## Use biomart's result to get the gene lengths etc.
dmel_annotations <- sm(get_biomart_annotations(species="dmelanogaster"))
## And ontology cateogies.
dmel_ontologies <- sm(get_biomart_ontologies(species="dmelanogaster"))
## Get the annotations ready to be recast as a gff file.
dmel_annotations$strand <- ifelse(dmel_annotations$strand == "1", "+", "-")
## Then make them into a granges object
dmel_granges <- GenomicRanges::makeGRangesFromDataFrame(dmel_annotations, keep.extra.columns=TRUE)
## I got a weird error when the column was Type and not type, I suspect though that this line is not needed.
dmel_granges$type <- dmel_annotations$Type
## Recast the data frame first as a List of GRanges
dmel <- as.data.frame(dmel_granges)
dmel$ID <- dmel$geneID
## Get the gene lengths for goseq.
dmel_lengths <- dmel_annotations[, c("geneID", "length")]
colnames(dmel_lengths) <- c("ID","width")
rownames(dmel_lengths) <- make.names(dmel_lengths[["ID"]], unique=TRUE)
## Drop all duplicate gene IDs
dmel_lengths <- dmel_lengths[ !grepl("\\.", rownames(dmel_lengths)), ]

limma <- new.env()
load("de_limma.rda", envir=limma)
table <- limma$hpgl_table
z_sig_genes <- sm(get_sig_genes(table, column="untreated", z=1)$up_genes)
fc_sig_genes <- sm(get_sig_genes(table, column="untreated", fc=1)$up_genes)
fcp_sig_genes <- sm(get_sig_genes(table, column="untreated", fc=1, p=0.05)$up_genes)
top200_sig_genes <- sm(get_sig_genes(table, column="untreated", n=200)$up_genes)

expected <- 1286
actual <- nrow(z_sig_genes)
test_that("Do we get a consistent number of z-scored up genes?", {
    expect_equal(expected, actual)
})

expected <- 129
actual <- nrow(fc_sig_genes)
test_that("Do we get a consistent number of fc up genes?", {
    expect_equal(expected, actual)
})

##expected <- 118
## The new version of limma gets 116 instead of 118 now.
expected <- 116
actual <- nrow(fcp_sig_genes)
test_that("Do we get a consistent number of fc and p up genes?", {
    expect_equal(expected, actual)
})

expected <- 200
actual <- nrow(top200_sig_genes)
test_that("Do we get a consistent top 200 up genes?", {
    expect_equal(expected, actual)
})

expected <- c("FBgn0085359", "FBgn0039155", "FBgn0024288",
              "FBgn0039827", "FBgn0034434", "FBgn0034736")
## This is the first 6 entries in all of the following tables.

actual <- head(rownames(z_sig_genes))
test_that("Do we get a consistent top set of z-scored up genes?", {
    expect_equal(expected, actual)
})

actual <- head(rownames(fc_sig_genes))
test_that("Do we get a consistent top set of fc up genes?", {
    expect_equal(expected, actual)
})

actual <- head(rownames(fcp_sig_genes))
test_that("Do we get a consistent top set of z-scored up genes?", {
    expect_equal(expected, actual)
})

actual <- head(rownames(top200_sig_genes))
test_that("Do we get a consistent top set of z-scored up genes?", {
    expect_equal(expected, actual)
})

save(list=ls(), file="gsea_siggenes.rda")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 51gsea_siggenes.R in ", elapsed,  " seconds."))
