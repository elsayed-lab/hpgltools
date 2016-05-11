library(testthat)
library(hpgltools)
context("Test usability of DESeq2")

## This section is copy/pasted to all of these tests, that is dumb.
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts),]
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single-end","single-end","paired-end",
        "paired-end","single-end","paired-end","paired-end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

message("Performing DESeq2 differential expression analysis as per the DESeq vignette.")
summarized <- DESeq2::DESeqDataSetFromMatrix(countData=counts,
                                             colData=metadata,
                                             design=~ condition + batch)
dataset <- suppressMessages(DESeq2::DESeqDataSet(se=summarized, design=~ condition + batch))
deseq_sf <- suppressMessages(DESeq2::estimateSizeFactors(dataset))
deseq_disp <- suppressMessages(DESeq2::estimateDispersions(deseq_sf))
deseq_run = suppressMessages(DESeq2::nbinomWaldTest(deseq_disp))
deseq_result <- as.data.frame(DESeq2::results(deseq_run,
                                              contrast=c("condition", "treated", "untreated"),
                                              format="DataFrame"))


message("Performing DESeq2 analysis using hpgltools.")
pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)
hpgl_deseq <- suppressMessages(deseq_pairwise(pasilla_expt, model_batch=TRUE))
## Note that running the all_pairwise family of functions results in arbitrarily chosen x/y which may be
## the opposite of what you actually want.
## Also, the all_pairwise functions reorders the result by logFC.
## In addition, I force it to a limited number of significant digits to avoid have data structures with 35 decimals.
hpgl_result <- hpgl_deseq[["all_tables"]][["untreated_vs_treated"]]
hpgl_result$logFC <- (-1 * hpgl_result$logFC)
hpgl_result_reordered <- hpgl_result[order(hpgl_result[["logFC"]]),]
deseq_result_reordered <- deseq_result[order(deseq_result[["log2FoldChange"]]),]

hpgl_logfc <- hpgl_result_reordered$logFC
deseq_logfc <- deseq_result_reordered$log2FoldChange
hpgl_basemean <- hpgl_result_reordered$baseMean
deseq_basemean <- deseq_result_reordered$baseMean
hpgl_stat <- hpgl_result_reordered$stat
deseq_stat <- deseq_result_reordered$stat
test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise()?", {
          expect_equal(deseq_logfc, hpgl_logfc, tolerance=0.001)
          expect_equal(deseq_basemean, hpgl_basemean, tolerance=1)
          expect_equal(deseq_stat, hpgl_stat * -1, tolerance=0.1)
})

message("Finished tests in 24deseq.")
