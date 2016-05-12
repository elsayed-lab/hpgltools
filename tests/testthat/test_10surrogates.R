library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)

context("Do surrogate estimators provide expected outputs?")

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts),]
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

counts <- counts[rowSums(counts) > ncol(counts),]

pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata)

pasilla_surrogates <- suppressMessages(compare_surrogate_estimates(pasilla_expt))
pca_char_adjust <- as.numeric(c("-0.108919503928491", "-0.570855719628606",
                     "0.413293224481248", "0.293317999343559",
                     "-0.527145619820086", "0.284623554215443",
                     "0.215686065336933"))
test_that("Does the compare_surrogate stuff work?", {
    expect_equal(pca_char_adjust,
                 as.numeric(as.character(pasilla_surrogates[["pca_adjust"]][["model_adjust"]])),
                 tolerance = 0.0001)
})
