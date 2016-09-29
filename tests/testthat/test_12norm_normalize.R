library(testthat)
library(hpgltools)
context("Are normalizations consistent over time? (Normalizations)")

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
    libType=c("single_end","single_end","paired_end",
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


## Test normalizations -- I should change this to be automatically generated for expected
expected <- as.numeric(c(5.857143, 91.500000, 4400.000000, 543.785714, 10.714286))
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="quant"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- as.numeric(actual_df[test_genes, c("untreated1")])
test_that("quant normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Similar test for size-factor normalization
expected <- c(4.392658, 80.824908, 4097.471407, 512.183926, 8.785316)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="sf"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("size-factor normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Check another size-factor normalization
expected <- c(4.392658, 80.824908, 4097.471407, 512.183926, 8.785316)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="sf2"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("size-factor2 normalization gives expected values?", {
    expect_equal(expected, actual)
})

## Oh I never noticed before that this is a log, too
expected <- c(5.488150, 7.082043, 12.021996, 9.160395, 5.707992)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="vsd"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("vsd normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## FIXME
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qsmooth")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink")
##pasilla_norm <- normalize_expt(pasilla_expt, norm="qshrink_median")

expected <- c(4.927997, 91.830657, 4765.366532, 613.466245, 9.342734)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="tmm"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("tmm normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.902692, 90.336774, 4803.090308, 608.726226, 9.488822)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="upperquartile"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("upperquartile normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(4.927854, 91.079703, 4840.296148, 615.582521, 9.205998)
names(expected) <- test_genes
pasilla_norm <- sm(normalize_expt(pasilla_expt, norm="rle"))
actual_df <- Biobase::exprs(pasilla_norm[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("RLE normalization gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

message("\n")
message("Finished 12norm_normalize.R")
