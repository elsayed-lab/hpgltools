library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Does pasilla load into hpgltools?")

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
metadata$sampleid <- rownames(metadata)

## Make sure it is still possible to create an expt
pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata, savefile="pasilla")
count_data <- as.matrix(counts)
hpgl_data <- Biobase::exprs(pasilla_expt$expressionset)
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(count_data, hpgl_data)
})

## Test that the expt has a design which makes sense.
known_samples <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
expt_samples <- as.character(pasilla_expt[["design"]][["sampleid"]])
known_conditions <- c("untreated","untreated","untreated","untreated","treated","treated","treated")
expt_conditions <- as.character(pasilla_expt[["design"]][["condition"]])
known_batches <- c("single_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end")
expt_batches <-  as.character(pasilla_expt[["design"]][["batch"]])
test_that("Is the experimental design maintained?", {
    expect_equal(known_samples, expt_samples)
    expect_equal(known_conditions, expt_conditions)
    expect_equal(known_batches, expt_batches)
})

known_libsizes <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(known_libsizes) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
expt_libsizes <- pasilla_expt[["libsize"]]
test_that("Are the library sizes intact?", {
    expect_equal(known_libsizes, expt_libsizes)
})
