start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
context("301load_data.R: Does pasilla load into hpgltools?\n")

## Try loading some annotation information for this species.

## This now generates an error on travis, but not on my computer.
gene_info <- sm(load_biomart_annotations(
    species="dmelanogaster", overwrite=TRUE))[["annotation"]]
info_idx <- gene_info[["gene_biotype"]] == "protein_coding"
gene_info <- gene_info[info_idx, ]
rownames(gene_info) <- make.names(gene_info[["ensembl_gene_id"]], unique=TRUE)

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts), ]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(
  row.names=colnames(counts),
  condition=c("untreated", "untreated", "untreated",
              "untreated", "treated", "treated", "treated"),
  libType=c("single_end", "single_end", "paired_end",
            "paired_end", "single_end", "paired_end", "paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
save(list=ls(), file="pasilla_df.rda")

## Make sure it is still possible to create an expt
pasilla_expt <- create_expt(count_dataframe=counts, metadata=metadata,
                            savefile="pasilla.rda", gene_info=gene_info)
## Recent changes to how my expressionsets are created mean that the order of
## genes is hard-set to the order of annotations in the annotation data and
## therefore _not_ the order of genes found in the count tables.
actual <- exprs(pasilla_expt)
actual <- actual[ order(row.names(actual)), ]
expected <- as.matrix(counts)
expected <- expected[ order(row.names(expected)), ]
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(expected, actual)
})

## The set of annotations should be in a consistent order.
hpgl_annotations <- fData(pasilla_expt)
chosen_genes <- c("FBgn0000008", "FBgn0000014", "FBgn0000017",
                  "FBgn0000018", "FBgn0000024", "FBgn0000032")

expected <- 10100
actual <- nrow(exprs(pasilla_expt))
test_that("Was the annotation information imported into the expressionset? (static rownames?)", {
    expect_gt(actual, expected)
})

## Then lengths of features should therefore remain consistent.
##expected <- c(1521, 192, 1344, 1428, 1428, 1428)
expected <- c(3987, 990, 4860, 1617, 1947, 1314)
actual <- as.numeric(hpgl_annotations[chosen_genes, "cds_length"])
##  head(sm(sort(as.numeric(hpgl_annotations[["cds_length"]]))))
test_that("Was the annotation information imported into the expressionset? (static lengths?)", {
    expect_equal(expected, actual)
})

## By the same token, the start positions of genes should remain consistent.
##expected <- c(18822604, 30212156, 7782797, 19116483, 19116483, 19116483)
expected <- c(22136968, 16807214, 16615866, 10973443, 13222951, 29991144)
actual <- as.numeric(hpgl_annotations[chosen_genes, "start_position"])
test_that("Was the annotation information imported into the expressionset? (static starts?)", {
    expect_equal(expected, actual)
})

## As should the chromosome arms of these genes.
##expected <- c("3L", "3R", "2R", "2L", "2L", "2L")
expected <- c("2R", "3R", "3L", "2L", "3R", "3R")
actual <- hpgl_annotations[chosen_genes, "chromosome_name"]
test_that("Was the annotation information imported into the expressionset? (static chromosomes?)", {
    expect_equal(expected, actual)
})

## Test that the expt has a design which makes sense.
expected <- c("untreated1","untreated2","untreated3","untreated4","treated1","treated2","treated3")
actual <- as.character(pasilla_expt[["design"]][["sampleid"]])
test_that("Is the experimental design maintained for samples?", {
    expect_equal(expected, actual)
})

## The conditions specified by the pasilla data set are treated and untreated and should not change.
expected <- c("untreated","untreated","untreated","untreated","treated","treated","treated")
actual <- as.character(pasilla_expt[["design"]][["condition"]])
test_that("Is the experimental design maintained for conditions?", {
    expect_equal(expected, actual)
})

## Some sequencing runs of pasilla are paired, and some are single ended; this should not change.
expected <- c("single_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end")
actual <-  as.character(pasilla_expt[["design"]][["batch"]])
test_that("Is the experimental design maintained for batches?", {
    expect_equal(expected, actual)
})

## We know a priori the library sizes, make sure that this information is intact.
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
names(expected) <- c("untreated1", "untreated2", "untreated3", "untreated4",
                     "treated1", "treated2", "treated3")
actual <- pasilla_expt[["libsize"]]
test_that("Are the library sizes intact?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 02load_data.R in ", elapsed,  " seconds."))
