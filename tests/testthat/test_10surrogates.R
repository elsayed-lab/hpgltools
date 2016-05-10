library(testthat)
library(hpgltools)
context("Testing surrogate testing.")
require.auto("pasilla")
library(pasilla)
data(pasillaGenes)

datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(
    "row.names" = colnames(counts),
    "condition" = c("untreated","untreated","untreated",
                    "untreated","treated","treated","treated"),
    "libType" = c("single-end","single-end","paired-end",
                  "paired-end","single-end","paired-end","paired-end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata$sampleid <- rownames(metadata)
## Make sure it is still possible to create an expt
message("Setting up an expt class to contain the pasilla data and metadata.")
pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata)

pasilla_surrogates <- compare_surrogate_estimates(pasilla_expt)

pca_char_adjust <- as.numeric(c("-0.108919503928491", "-0.570855719628606",
                     "0.413293224481248", "0.293317999343559",
                     "-0.527145619820086", "0.284623554215443",
                     "0.215686065336933"))
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(pca_char_adjust,
                 as.numeric(as.character(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])),
                 tolerance = 0.0001)
})

message("Hey, write more tests here.")
message("YAY! Finished testing surrogate estimators!")
