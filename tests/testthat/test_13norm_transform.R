library(testthat)
library(hpgltools)
context("Are normalizations consistent over time (Tranformations)?")

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

## Test transformations
expected <- c(2.584963, 6.539159, 12.187661, 9.189825, 3.459432)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform="log2"))
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("log2 transformation gives expected values?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(0.7781513, 1.9684829, 3.6688516, 2.7664128, 1.0413927)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform="log10"))
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("log10 transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(1.791759, 4.532599, 8.447843, 6.369901, 2.397895)
names(expected) <- test_genes
pasilla_trans <- sm(normalize_expt(pasilla_expt, transform="log"))
actual_df <- Biobase::exprs(pasilla_trans[["expressionset"]])
actual <- actual_df[test_genes, c("untreated1")]
test_that("loge transformation gives expected values (why log10!?)?", {
    expect_equal(expected, actual, tolerance=0.0001)
})
