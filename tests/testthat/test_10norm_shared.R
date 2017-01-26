start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("10norm_shared.R: Are normalizations consistent over time (Shared functions)?\n")

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

## First make sure the pasilla_expt still has the stuff we expect
expected <- "This is an expt class."
actual <- pasilla_expt[["title"]]
test_that("Pasilla title?", {
    expect_equal(expected, actual)
})

## Ensure that the beginning count table library sizes are identical.
expected <- colSums(counts)
actual <- pasilla_expt[["original_libsize"]]
names(expected) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
test_that("Pasilla libsize?", {
    expect_equal(expected, actual)
})

## Check a few arbitrary counts to make sure they are maintained.
expected <- counts[test_genes, "untreated1"]
testing_counts <- Biobase::exprs(pasilla_expt[["expressionset"]])
actual <- as.numeric(testing_counts[test_genes, "untreated1"])
test_that("Pasilla count tables? (untreated1)", {
    expect_equal(expected, actual)
})

## Check that all samples agree for 1 gene.
test_gene <- "FBgn0062565"
expected <- as.numeric(counts[test_gene, ])
actual <- as.numeric(Biobase::exprs(pasilla_expt[["expressionset"]])[test_gene, ])
expected <- c(4, 7, 3, 3, 9, 10, 9)
test_that("Pasilla count tables? (gene FBgn0063565)", {
    expect_equal(expected, actual)
})

## Ensure that normalize_expt does not mess up the data when called without arguments (this wasn't true once)
unmolested <- sm(normalize_expt(pasilla_expt))
expected <- as.matrix(Biobase::exprs(pasilla_expt[["expressionset"]]))  ## I accidently changed this to potentially return a data.frame
actual <- Biobase::exprs(unmolested[["expressionset"]])
test_that("Pasilla (un)normalized counts?", {
    expect_equal(expected, actual)
})

## Make sure that the pData information is maintained through normalization
expected <- Biobase::pData(pasilla_expt[["expressionset"]])
actual <- Biobase::pData(unmolested[["expressionset"]])
test_that("Pasilla (un)normalized pdata?", {
    expect_equal(expected, actual)
})

## Also ensure that the library sizes (which are very important for limma) are not messed up.
expected <- pasilla_expt[["libsize"]]
actual <- unmolested[["libsize"]]
test_that("Pasilla (un)normalized libsize?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 10norm_shared.R in ", elapsed, "seconds."))
