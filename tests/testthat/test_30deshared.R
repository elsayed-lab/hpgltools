library(testthat)
library(hpgltools)

context("Do the combined differential expression searches work?")

deseq <- new.env()
edger <- new.env()
limma <- new.env()
basic <- new.env()

load("de_deseq.rda", envir=deseq)
load("de_edger.rda", envir=edger)
load("de_limma.rda", envir=limma)
load("de_basic.rda", envir=basic)

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts), ]
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata$Sample.id <- rownames(metadata)

pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata)
normalized_expt <- normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm")

hpgl_result <- suppressMessages(all_pairwise(normalized_expt, model_batch=TRUE))

previous_deseq <- deseq$hpgl_deseq$all_tables[["untreated_vs_treated"]]
previous_edger <- edger$hpgl_edger$all_tables[["untreated_vs_treated"]]
previous_limma <- limma$hpgl_limma$all_tables[["untreated_vs_treated"]]
previous_basic <- basic$hpgl_basic$all_tables[["untreated_vs_treated"]]

this_deseq <- hpgl_result$deseq$all_tables$untreated_vs_treated
this_edger <- hpgl_result$edger$all_tables$untreated_vs_treated
this_limma <- hpgl_result$limma$all_tables$untreated_vs_treated
this_basic <- hpgl_result$basic$all_tables$untreated_vs_treated

d <- e <- b <- l <- NULL
d <- expect_equal(previous_deseq, this_deseq)
e <- expect_equal(previous_edger, this_edger)
l <- expect_equal(previous_limma, this_limma)
b <- expect_equal(previous_basic, this_basic)
test_that("Do we get similar results to previous DE runs?", {
    expect_equal(previous_deseq, this_deseq)
    expect_equal(previous_edger, this_edger)
    expect_equal(previous_limma, this_limma)
    expect_equal(previous_basic, this_basic)
})

if (is.null(d) | is.null(e) | is.null(l) | is.null(b)) {
    skip("Debugging these output files, not deleting the .rda files.")
} else {
    le <- hpgl_result$comparison$comp[[1]]
    ld <- hpgl_result$comparison$comp[[2]]
    ed <- hpgl_result$comparison$comp[[3]]
    lb <- hpgl_result$comparison$comp[[4]]
    eb <- hpgl_result$comparison$comp[[5]]
    db <- hpgl_result$comparison$comp[[6]]
    test_that("Are the comparisons between DE tools sufficiently similar?", {
        expect_gt(le, 0.95)
        expect_gt(ld, 0.80)
        expect_gt(ed, 0.75)
        expect_gt(lb, 0.88)
        expect_gt(eb, 0.83)
        expect_gt(db, 0.68)
    })
    de_removed <- file.remove("de_deseq.rda")
    ed_removed <- file.remove("de_edger.rda")
    li_removed <- file.remove("de_limma.rda")
    ba_removed <- file.remove("de_basic.rda")
}
