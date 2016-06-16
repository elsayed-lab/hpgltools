library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Does pasilla load into hpgltools?")

## Try loading some annotation information for this species.
gene_info <- s_p(get_biomart_annotations(species="dmelanogaster"))[["result"]]
info_idx <- gene_info[["Type"]] == "protein_coding"
gene_info <- gene_info[info_idx, ]
rownames(gene_info) <- make.names(gene_info[["geneID"]], unique=TRUE)

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

## Make sure it is still possible to create an expt
pasilla_expt <- create_expt(count_dataframe=counts, metadata=metadata, savefile="pasilla", gene_info=gene_info)
## Recent changes to how my expressionsets are created mean that the order of genes is hard-set to the order of annotations
## in the annotation data and therefore _not_ the order of genes found in the count tables.
actual <- as.matrix(Biobase::exprs(pasilla_expt[["expressionset"]]))
actual <- actual[ order(row.names(actual)), ]
expected <- as.matrix(counts)
expected <- expected[ order(row.names(expected)), ]
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(expected, actual)
})

hpgl_annotations <- Biobase::fData(pasilla_expt[["expressionset"]])
actual <- head(rownames(hpgl_annotations))
expected <- c("FBgn0260439", "FBgn0031081", "FBgn0062565", "FBgn0031089", "FBgn0031092", "FBgn0031094")
test_that("Was the annotation information imported into the expressionset? (static rownames?)", {
    expect_equal(expected, actual)
})
actual <- head(hpgl_annotations[["length"]])
expected <- c(1776, 2361, 1164, 1326, 939, 819)
test_that("Was the annotation information imported into the expressionset? (static lengths?)", {
    expect_equal(expected, actual)
})
actual <- head(hpgl_annotations[["start"]])
expected <- c(8366038, 19961297, 20094398, 20148124, 20166673, 20170222)
test_that("Was the annotation information imported into the expressionset? (static starts?)", {
    expect_equal(expected, actual)
})
actual <- head(hpgl_annotations[["chromosome"]])
expected <- c("2L", "X", "X", "X", "X", "X")
test_that("Was the annotation information imported into the expressionset? (static chromosomes?)", {
    expect_equal(expected, actual)
})

## Test that the expt has a design which makes sense.
actual <- as.character(pasilla_expt[["design"]][["sampleid"]])
expected <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
test_that("Is the experimental design maintained for samples?", {
    expect_equal(expected, actual)
})

actual <- as.character(pasilla_expt[["design"]][["condition"]])
expected <- c("untreated","untreated","untreated","untreated","treated","treated","treated")
test_that("Is the experimental design maintained for conditions?", {
    expect_equal(expected, actual)
})

actual <-  as.character(pasilla_expt[["design"]][["batch"]])
expected <- c("single_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end")
test_that("Is the experimental design maintained for batches?", {
    expect_equal(expected, actual)
})

actual <- pasilla_expt[["libsize"]]
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(expected) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
test_that("Are the library sizes intact?", {
    expect_equal(expected, actual)
})
