library(testthat) 
library(hpgltools)
context("Test usability of EdgeR")
autoloads_all()
require.auto("pasilla")
require.auto("edgeR")

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

raw = DGEList(counts=counts, group=metadata$condition)
norm = calcNormFactors(raw)
disp_norm = estimateCommonDisp(norm)
tagdispnorm = estimateTagwiseDisp(disp_norm)
exact_test = exactTest(tagdispnorm)
exact_result = as.data.frame(topTags(exact_test, n=nrow(raw)))

exact_vs_cbcb = merge(cbcb_top, exact_result, by.x="row.names", by.y="row.names")
cor.test(exact_vs_cbcb$logFC.x, exact_vs_cbcb$logFC.y)

condition = design$condition
edger_design = model.matrix(~0 + condition)

glm_norm = estimateGLMCommonDisp(tagdispnorm, edger_design)
glm_trended = estimateGLMTrendedDisp(glm_norm, edger_design)
glm_tagged = estimateGLMTagwiseDisp(glm_trended, edger_design)
edger_fit = glmFit(glm_tagged, design=edger_design)
lrt = glmLRT(edger_fit, coef=2)
glm_result = as.data.frame(topTags(lrt, n=nrow(raw)))

glm_vs_exact = merge(exact_result, glm_result, by.x="row.names", by.y="row.names")
cor.test(glm_vs_exact$logFC.x, glm_vs_exact$logFC.y)

hpgl_edger = edger_pairwise(expt=pasilla_expt)
hpglglm_result = hpgl_edger$results$untreated_minus_treated
hpgl_vs_exact = merge(exact_result, hpglglm_result, by.x="row.names", by.y="row.names")
cor.test(hpgl_vs_exact$logFC.x, glm_vs_exact$logFC.y)
