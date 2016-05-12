library(testthat)
library(hpgltools)

context("Does the basic differential expression analysis work?")

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata$Sample.id <- rownames(metadata)

pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata)
norm_expt <- suppressMessages(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm"))

hpgl_pas_basic <- suppressMessages(basic_pairwise(pasilla_expt))
hpgl_basic <- suppressMessages(basic_pairwise(norm_expt))

test_that("Does a non-normalized basic run equal a normalized basic run?", {
    expect_equal(hpgl_pas_basic$all_tables$untreated_vs_treated, hpgl_basic$all_tables$untreated_vs_treated)
})

expected_medians_treated <- c(2.8356541, -4.0744928, 8.3510620, 5.0276022, -0.9340304, 6.2167638)
actual_medians_treated <- head(hpgl_basic[["medians"]][["treated"]])
expected_medians_untreated <- c(2.919061, -4.074470, 8.513323, 5.099943, -1.109891, 6.266102)
actual_medians_untreated <- head(hpgl_basic[["medians"]][["untreated"]])
expected_logfc <- c(8.341e-02, 2.234e-05, 1.623e-01, 7.234e-02, -1.759e-01, 4.934e-02)
actual_logfc <- head(hpgl_basic[["all_tables"]][["untreated_vs_treated"]][["logFC"]])

test_that("Do we get the values for treated and untreated samples?", {
    expect_equal(expected_medians_treated, actual_medians_treated, tolerance=0.001)
    expect_equal(expected_medians_untreated, actual_medians_untreated, tolerance=0.001)
    expect_equal(expected_logfc, actual_logfc, tolerance=0.001)
})

save(list=ls(), file="de_basic.rda")
