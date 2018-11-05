start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("24de_deseq.R: Does hpgltools work with DESeq2?\n")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("de_limma.rda", envir=limma)
counts <- limma[["counts"]]
design <- limma[["design"]]

metadata <- design
colnames(metadata) <- c("condition", "batch")
## Performing DESeq2 differential expression analysis as per the DESeq vignette.
summarized <- DESeq2::DESeqDataSetFromMatrix(countData=counts,
                                             colData=metadata,
                                             design=~ 0 + batch + condition)
dataset <- DESeq2::DESeqDataSet(se=summarized, design=~ 0 + batch + condition)
deseq_sf <- DESeq2::estimateSizeFactors(dataset)
deseq_disp <- DESeq2::estimateDispersions(deseq_sf)
deseq_run <- DESeq2::nbinomWaldTest(deseq_disp, betaPrior=FALSE)
deseq_result <- as.data.frame(DESeq2::results(deseq_run,
                                              contrast=c("condition", "treated", "untreated"),
                                              format="DataFrame"))

## Performing DESeq2 analysis using hpgltools.
hpgl_deseq <- sm(deseq2_pairwise(input=pasilla_expt,
                                 model_batch=TRUE,
                                 deseq_excel="deseq_test.xlsx"))
test_that("Can I write a deseq2 table?", {
    expect_true(file.exists("deseq_test.xlsx"))
})

## Note that running the all_pairwise family of functions results in arbitrarily
## chosen x/y which may be the opposite of what you actually want.
## Also, the all_pairwise functions reorders the result by logFC.
## In addition, I force it to a limited number of significant digits to avoid
## have data structures with 35 decimals. This is explicitly in the opposite
## order as what deseq_pairwise will choose in order to make it necessary to
## test the columns that change as a result.
hpgl_result <- hpgl_deseq[["all_tables"]][["untreated_vs_treated"]]
hpgl_result_reordered <- hpgl_result[order(rownames(hpgl_result)), ]
deseq_result_reordered <- deseq_result[order(rownames(deseq_result)), ]

## Columns to test: baseMean, logFC, lfcSE, stat, P.Value, adj.P.Val, qvalue
deseq_basemean <- deseq_result_reordered[["baseMean"]]
hpgl_basemean <- hpgl_result_reordered[["baseMean"]]
deseq_logfc <- deseq_result_reordered[["log2FoldChange"]]
hpgl_logfc <- hpgl_result_reordered[["logFC"]] * -1
deseq_lfcse <- deseq_result_reordered[["lfcSE"]]
hpgl_lfcse <- hpgl_result_reordered[["lfcSE"]]
deseq_stat <- deseq_result_reordered[["stat"]]
hpgl_stat <- hpgl_result_reordered[["stat"]] * -1
deseq_pval <- deseq_result_reordered[["pvalue"]]
hpgl_pval <- hpgl_result_reordered[["P.Value"]]
deseq_adjpval <- deseq_result_reordered[["adj.P.Val"]]
hpgl_adjpval <- hpgl_result_reordered[["padj"]]

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): basemean?", {
    expect_equal(deseq_basemean, hpgl_basemean, tolerance=1)
})

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): logfc?", {
    expect_equal(deseq_logfc, hpgl_logfc, tolerance=0.001)
})

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): lfcse?", {
    expect_equal(deseq_lfcse, hpgl_lfcse, tolerance=0.2)
})

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): stat?", {
    expect_equal(deseq_stat, hpgl_stat, tolerance=0.1)
})

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): pval?", {
    expect_equal(deseq_pval, hpgl_pval, tolerance=0.1)
})

test_that("Does the DESeq2 vignette agree with the result from deseq_pairwise(): adjpval?", {
    expect_equal(deseq_adjpval, hpgl_adjpval, tolerance=0.1)
})

save(list=ls(), file="de_deseq.rda")
end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 24de_deseq.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
