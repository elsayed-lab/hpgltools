library(testthat)
library(hpgltools)

context("Does edgeR work with hpgltools?")

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
pair <- "untreated - treated"
contr <- makeContrasts(contrasts=pair, levels=model)
glm_result <- edgeR::glmLRT(glmfit, contrast=contr)
glm_table <- as.data.frame(edgeR::topTags(glm_result, n=nrow(raw), sort.by="logFC"))

## Create the expt object
pasilla_expt <- create_expt(count_dataframe=counts, meta_dataframe=metadata)
cbcb_data <- as.matrix(counts)
hpgl_data <- Biobase::exprs(pasilla_expt$expressionset)
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(cbcb_data, hpgl_data)
})

## Perform the edgeR analysis in hpgltools
hpgl_edger <- suppressMessages(edger_pairwise(pasilla_expt))
hpgl_result <- hpgl_edger$all_tables$untreated_vs_treated
hpgl_reordered <- hpgl_result[order(hpgl_result[["logFC"]]), ]
edger_reordered <- glm_table[order(glm_table[["logFC"]]), ]

## Compare the edgeR analysis to that in hpgltools
test_that("Is the hpgl pairwise similar to edgeR's default method?", {
    expect_equal(hpgl_reordered$logFC, edger_reordered$logFC, tolerance=0.1)
    expect_equal(hpgl_reordered$PValue, edger_reordered$PValue, tolerance=0.1)
})

save(list=ls(), file="de_edger.rda")
