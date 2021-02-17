start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("351gsea_siggenes.R: Do we get consistent sets of 'significant' genes for GSEA analyses?
  1234567890\n")

limma <- new.env()
load("320_de_limma.rda", envir = limma)
table <- limma[["hpgl_limma"]][["all_tables"]][[1]]

## Use biomart's result to get the gene lengths etc.
dmel_annotations <- sm(load_biomart_annotations(species = "dmelanogaster"))
dmel_annotations <- dmel_annotations[["annotation"]]
## And ontology cateogies.
dmel_ontologies <- sm(load_biomart_go(species = "dmelanogaster"))
dmel_ontologies <- dmel_ontologies[["go"]]
## Get the annotations ready to be recast as a gff file.
dmel_annotations[["strand"]] <- ifelse(dmel_annotations[["strand"]] == "1", "+", "-")
colnames(dmel_annotations) <- c("transcript_id", "gene_id", "description", "type",
                                "length", "chromosome", "strand",
                                "start", "end")
## Then make them into a granges object
## WTF is with 'unable to find an inherited method for function 'seqinfo' for signature "GRanges". ??
## Why is it that if I start with a fresh R session, this works fine!?!?!?

dmel_granges <- GenomicRanges::makeGRangesFromDataFrame(dmel_annotations, keep.extra.columns = TRUE)
## I got a weird error when the column was Type and not type, I suspect though that this line is not needed.
dmel_granges$Type <- dmel_annotations[["type"]]
## Recast the data frame first as a List of GRanges
dmel <- as.data.frame(dmel_granges)
dmel[["ID"]] <- dmel[["gene_id"]]
## Get the gene lengths for goseq.
dmel_lengths <- dmel_annotations[, c("gene_id", "length")]
colnames(dmel_lengths) <- c("ID","width")
rownames(dmel_lengths) <- make.names(dmel_lengths[["ID"]], unique = TRUE)
## Drop all duplicate gene IDs
dmel_lengths <- dmel_lengths[ !grepl("\\.", rownames(dmel_lengths)), ]

z_sig_genes <- sm(get_sig_genes(table, column = "logFC", z = 1)[["up_genes"]])
fc_sig_genes <- sm(get_sig_genes(table, column = "logFC", lfc = 1)[["up_genes"]])
fcp_sig_genes <- sm(get_sig_genes(table, column = "logFC", lfc = 1, p = 0.05)[["up_genes"]])
top200_sig_genes <- sm(get_sig_genes(table, column = "logFC", n = 200)[["up_genes"]])

expected <- 1852
actual <- nrow(z_sig_genes)
test_that("Do we get a consistent number of z-scored up genes?", {
    expect_equal(expected, actual)
})

expected <- 331
actual <- nrow(fc_sig_genes)
test_that("Do we get a consistent number of fc up genes?", {
    expect_equal(expected, actual)
})

##expected <- 118
## The new version of limma gets 116 instead of 118 now.
expected <- 123
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

save(list = ls(), file = "351_gsea_siggenes.rda")
end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 351gsea_siggenes.R in ", elapsed,  " seconds.")
