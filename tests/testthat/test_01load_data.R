start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("01load_data.R: Does pasilla load into hpgltools?\n")

## Try loading some annotation information for this species.
gene_info <- sm(get_biomart_annotations(species="dmelanogaster"))
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
pasilla_expt <- sm(create_expt(count_dataframe=counts, metadata=metadata, savefile="pasilla", gene_info=gene_info))
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
expected <- c("FBgn0000008", "FBgn0000014", "FBgn0000017", "FBgn0000018", "FBgn0000024", "FBgn0000032")
actual <- head(sort(rownames(hpgl_annotations)))
test_that("Was the annotation information imported into the expressionset? (static rownames?)", {
    expect_equal(expected, actual)
})

expected <- c(78, 81, 99, 123, 123, 123)
actual <- head(sort(hpgl_annotations[["length"]]))
test_that("Was the annotation information imported into the expressionset? (static lengths?)", {
    expect_equal(expected, actual)
})

expected <- c(7529, 9839, 21823, 25402, 32478, 47710)
actual <- head(sort(hpgl_annotations[["start"]]))
test_that("Was the annotation information imported into the expressionset? (static starts?)", {
    expect_equal(expected, actual)
})

expected <- c("2R", "3R", "3L", "2L", "3R", "3R")
actual <- head(hpgl_annotations[["chromosome"]])
test_that("Was the annotation information imported into the expressionset? (static chromosomes?)", {
    expect_equal(expected, actual)
})

## Test that the expt has a design which makes sense.
expected <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
actual <- as.character(pasilla_expt[["design"]][["sampleid"]])
test_that("Is the experimental design maintained for samples?", {
    expect_equal(expected, actual)
})

expected <- c("untreated","untreated","untreated","untreated","treated","treated","treated")
actual <- as.character(pasilla_expt[["design"]][["condition"]])
test_that("Is the experimental design maintained for conditions?", {
    expect_equal(expected, actual)
})

expected <- c("single_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end")
actual <-  as.character(pasilla_expt[["design"]][["batch"]])
test_that("Is the experimental design maintained for batches?", {
    expect_equal(expected, actual)
})

expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(expected) <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
actual <- pasilla_expt[["libsize"]]
test_that("Are the library sizes intact?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end - start), digits=1)
message(paste0("\nFinished 01load_data.R in ", elapsed,  " seconds."))
