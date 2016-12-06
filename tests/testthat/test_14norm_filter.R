start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("14norm_filter.R: Are normalizations consistent over time (Filtering)?\n")

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




## Test filter
expected <- c(7526, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter="cbcb"))
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("cbcb filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter="pofa"))
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("pofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter="kofa"))
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("kofa filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

## TODO These may need adjustment
expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter="cv"))
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("cv filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

expected <- c(10153, 7)
pasilla_filter <- sm(normalize_expt(pasilla_expt, filter="simple"))
actual <- dim(Biobase::exprs(pasilla_filter[["expressionset"]]))
test_that("simple filtering leaves behind the expected number of genes?", {
    expect_equal(expected, actual, tolerance=0.0001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 14norm_filter.R in ", elapsed,  " seconds."))
