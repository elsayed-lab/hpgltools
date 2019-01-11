start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
## The following library() seems to be required for glmQLFTest() to work.
tt <- library(edgeR)
context("26de_edger.R: Does hpgltools work with edgeR?\n")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("de_limma.rda", envir=limma)
counts <- limma[["counts"]]
design <- limma[["design"]]

metadata <- design
colnames(metadata) <- c("condition", "batch")
## Performing edgeR differential expression analysis as per the edgeR vignette.
model <- model.matrix(~ 0 + design[["condition"]] + design[["libType"]])
colnames(model) <- c("treated", "untreated", "libtype")
raw <- edgeR::DGEList(counts=counts, group=metadata[["condition"]])
norm <- edgeR::calcNormFactors(raw)
disp_norm <- edgeR::estimateCommonDisp(norm)
tagdispnorm <- edgeR::estimateTagwiseDisp(disp_norm)
glmnorm <- edgeR::estimateGLMCommonDisp(tagdispnorm, model)
glmtrend <- edgeR::estimateGLMTrendedDisp(glmnorm, model)
glmtagged <- edgeR::estimateGLMTagwiseDisp(glmtrend, model)
glmfit <- edgeR::glmQLFit(glmtagged, design=model, robust=TRUE)
## This is explicitly in the opposite order as what edger_pairwise will choose in order to
## make it necessary to test the columns that change as a result.
pair <- "treated - untreated"
contr <- limma::makeContrasts(contrasts=pair, levels=model)
## Why does this keep failing!?
glm_result <- edgeR::glmQLFTest(glmfit, contrast=contr)
glm_table <- as.data.frame(edgeR::topTags(glm_result, n=nrow(raw), sort.by="logFC"))

## Create the expt object
expected <- as.matrix(counts)
expected <- expected[sort(rownames(expected)), ]
actual <- exprs(pasilla_expt)
actual <- actual[sort(rownames(actual)), ]
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(expected, actual)
})

## Perform the edgeR analysis in hpgltools
hpgl_edger <- sm(edger_pairwise(pasilla_expt, edger_method="long", edger_test="qlf"))

hpgl_result <- hpgl_edger[["all_tables"]][["untreated_vs_treated"]]
hpgl_result[["logFC"]] <- hpgl_result[["logFC"]] * -1

## Because of rounding errors, the order of logCPM with respect to logFC is not
## maintained from hpgl->edger
## Therefore, order by rownames!
edger_reordered <- glm_table[order(rownames(glm_table)), ]
hpgl_reordered <- hpgl_result[order(rownames(hpgl_result)), ]

## Columns to check: logFC,logCPM,LR,PValue,FDR,qvalue vs. logFC,logCPM,LR,PValue,FDR
edger_logfc <- edger_reordered[["logFC"]]
hpgl_logfc <- hpgl_reordered[["logFC"]]
edger_logcpm <- edger_reordered[["logCPM"]]
hpgl_logcpm <- hpgl_reordered[["logCPM"]]
edger_f <- edger_reordered[["F"]]
hpgl_f <- hpgl_reordered[["F"]]
edger_pval <- edger_reordered[["PValue"]]
hpgl_pval <- hpgl_reordered[["PValue"]]
edger_fdr <- edger_reordered[["FDR"]]
hpgl_fdr <- hpgl_reordered[["FDR"]]

test_that("Is the hpgl pairwise similar to edgeR's default method (logfc)?", {
    expect_equal(edger_logfc, hpgl_logfc, tolerance=0.0001)
})

test_that("Is the hpgl pairwise similar to edgeR's default method (logcpm)?", {
    expect_equal(edger_logcpm, hpgl_logcpm, tolerance=0.0001)
})

test_that("Is the hpgl pairwise similar to edgeR's default method (F)?", {
    expect_equal(edger_f, hpgl_f, tolerance=0.0001)
})

test_that("Is the hpgl pairwise similar to edgeR's default method (pval)?", {
    expect_equal(edger_pval, hpgl_pval, tolerance=0.0001)
})

test_that("Is the hpgl pairwise similar to edgeR's default method (fdr)?", {
    expect_equal(edger_fdr, hpgl_fdr, tolerance=0.0001)
})

edger_written <- write_edger(hpgl_edger, excel="edger_test.xlsx")
test_that("Can we write the results of an edger pairwise analysis?", {
    expect_true(file.exists("edger_test.xlsx"))
})

hpgl_edger <- sm(edger_pairwise(pasilla_expt))
save(list=ls(), file="de_edger.rda")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 26de_edger.R in ", elapsed,  " seconds."))
tt <- try(clear_session())
