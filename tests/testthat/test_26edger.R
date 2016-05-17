library(testthat)
library(hpgltools)

context("Does edgeR work with hpgltools?")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("de_limma.rda", envir=limma)
counts <- limma$counts
design <- limma$design

metadata <- design
colnames(metadata) <- c("condition", "batch")
## Performing edgeR differential expression analysis as per the edgeR vignette.
model <- model.matrix(~ 0 + design$condition + design$libType)
colnames(model) <- c("treated","untreated","libtype")
raw <- edgeR::DGEList(counts=counts, group=metadata$condition)
norm <- edgeR::calcNormFactors(raw)
disp_norm <- edgeR::estimateCommonDisp(norm)
tagdispnorm <- edgeR::estimateTagwiseDisp(disp_norm)
glmnorm <- edgeR::estimateGLMCommonDisp(tagdispnorm, model)
glmtrend <- edgeR::estimateGLMTrendedDisp(glmnorm, model)
glmtagged <- edgeR::estimateGLMTagwiseDisp(glmtrend, model)
glmfit <- edgeR::glmFit(glmtagged, design=model)
## This is explicitly in the opposite order as what edger_pairwise will choose in order to
## make it necessary to test the columns that change as a result.
pair <- "treated - untreated"
contr <- limma::makeContrasts(contrasts=pair, levels=model)
glm_result <- edgeR::glmLRT(glmfit, contrast=contr)
glm_table <- as.data.frame(edgeR::topTags(glm_result, n=nrow(raw), sort.by="logFC"))

## Create the expt object
cbcb_data <- as.matrix(counts)
hpgl_data <- Biobase::exprs(pasilla_expt$expressionset)
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(cbcb_data, hpgl_data)
})

## Perform the edgeR analysis in hpgltools
hpgl_edger <- sp(edger_pairwise(pasilla_expt))$result

hpgl_result <- hpgl_edger$all_tables$untreated_vs_treated
hpgl_result[["logFC"]] <- hpgl_result[["logFC"]] * -1

## Because of rounding errors, the order of logCPM with respect to logFC is not maintained from hpgl->edger
## Therefore, order by rownames!
edger_reordered <- glm_table[order(rownames(glm_table)), ]
hpgl_reordered <- hpgl_result[order(rownames(hpgl_result)), ]

## Columns to check: logFC,logCPM,LR,PValue,FDR,qvalue vs. logFC,logCPM,LR,PValue,FDR
edger_logfc <- edger_reordered$logFC
hpgl_logfc <- hpgl_reordered$logFC
edger_logcpm <- edger_reordered$logCPM
hpgl_logcpm <- hpgl_reordered$logCPM
edger_lr <- edger_reordered$LR
hpgl_lr <- hpgl_reordered$LR
edger_pval <- edger_reordered$PValue
hpgl_pval <- hpgl_reordered$PValue
edger_fdr <- edger_reordered$FDR
hpgl_fdr <- hpgl_reordered$FDR

test_that("Is the hpgl pairwise similar to edgeR's default method?", {
    expect_equal(edger_logfc, hpgl_logfc, tolerance=0.1)
    expect_equal(edger_logcpm, hpgl_logcpm, tolerance=0.1)
    expect_equal(edger_lr, hpgl_lr, tolerance=0.1)
    expect_equal(edger_pval, hpgl_pval, tolerance=0.1)
    expect_equal(edger_fdr, hpgl_fdr, tolerance=0.1)
})

save(list=ls(), file="de_edger.rda")
