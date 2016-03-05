library(testthat)
library(hpgltools)
context("Test usability of EdgeR")
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

raw = edgeR::DGEList(counts=counts, group=metadata$condition)
norm = edgeR::calcNormFactors(raw)
disp_norm = edgeR::estimateCommonDisp(norm)
tagdispnorm = edgeR::estimateTagwiseDisp(disp_norm)
exact_test = edgeR::exactTest(tagdispnorm)
exact_result = as.data.frame(edgeR::topTags(exact_test, n=nrow(raw)))

## Load the pasilla data set
message("Loading pasilla, setting up count tables.")
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single-end","single-end","paired-end",
        "paired-end","single-end","paired-end","paired-end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

## Make sure it is still possible to create an expt
message("Setting up an expt class to contain the pasilla data and metadata.")
pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)
cbcb_data = counts
pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)
hpgl_data = exprs(pasilla_expt$expressionset)
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(cbcb_data, hpgl_data)
})

hpgl_edger = edger_pairwise(pasilla_expt)
hpgl_result = hpgl_edger$all_tables$untreated_vs_treated
hpgl_vs_edger = merge(exact_result, hpgl_result, by.x="row.names", by.y="row.names")
similarity <- as.numeric(cor.test(hpgl_vs_edger$logFC.x, hpgl_vs_edger$logFC.y)[[4]])
test_that("Is the hpgl pairwise similar to edgeR's default method?", {
    expect_gt(similarity, 0.999)
})


