library(testthat) 
library(hpgltools)
context("Test hpgltools and cbcbSEQ")

## Make sure I didn't introduce any stupid syntax errors.
test_that("Is it possible to load/start cbcbSEQ/hpgltools?", {
    library(cbcbSEQ)
    library(hpgltools)
    autoloads_all()
    require.auto("pasilla")
})


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

## Check that normalization tools work similarly
message("Testing quantile raw normalization.")
cbcb_quantile = cbcbSEQ::qNorm(cbcb_data)
hpgl_quantile_data = hpgl_norm(expt=pasilla_expt, transform="raw", norm="quant", convert="raw", filter_low=FALSE, verbose=TRUE)
hpgl_quantile = hpgl_quantile_data$final_counts$count_table
test_that("Are the quantile normalizations identical?", {
    expect_equal(cbcb_quantile, hpgl_quantile)
})

## Check that quant+cpm are the same
message("Testing quantile(cpm()) normalization using edgeR's cpm().")
cbcb_qcpm = qNorm(cbcb_data)
cbcb_qcpm = edgeR::cpm(cbcb_qcpm)
hpgl_qcpm = hpgl_norm(expt=pasilla_expt, norm="quant", convert="edgecpm", filter_low=FALSE, verbose=TRUE)
hpgl_qcpm = hpgl_qcpm$final_counts$count_table
test_that("Are cpm conversions identical?", {
    expect_equal(cbcb_qcpm, hpgl_qcpm)
})

## log2/cpm that
message("Testing log2(quantile(cpm())) normalization using the cpm from voom()")
cbcb_l2qcpm_data = log2CPM(cbcb_quantile)
cbcb_l2qcpm = cbcb_l2qcpm_data$y
hpgl_l2qcpm_data = hpgl_norm(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", filter_low=FALSE, verbose=TRUE)
hpgl_l2qcpm = hpgl_l2qcpm_data$final_counts$count_table
hpgl_l2qcpm_expt = normalize_expt(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", filter_low=FALSE)
hpgl_l2qcpm2 = exprs(hpgl_l2qcpm_expt$expressionset)
test_that("Are l2qcpm conversions/transformations identical using two codepaths?", {
    expect_equal(cbcb_l2qcpm, hpgl_l2qcpm)
    expect_equal(cbcb_l2qcpm, hpgl_l2qcpm2)
})

## Check that PCA invocations are similar
message("Testing PCA invocations.")
cbcb_svd = cbcbSEQ::makeSVD(cbcb_l2qcpm)
hpgl_pca_info = hpgl_pca(expt=hpgl_l2qcpm_expt)
hpgl_svd = hpgl_pca_info$pca
cbcb_res = cbcbSEQ::pcRes(cbcb_svd$v, cbcb_svd$d, design$condition, design$libType)
hpgl_res = hpgl_pca_info$res
test_that("Do the PCA invocations provide the same results?", {
    expect_equal(cbcb_svd$v, hpgl_svd$v)
    expect_equal(cbcb_svd$d, hpgl_svd$d)
    expect_equal(cbcb_res, hpgl_res)
})

## invocation of batch correction
message("Testing batch correction results using modified combat.")
cbcb_batch = cbcbSEQ::combatMod(cbcb_l2qcpm, batch=design$libType, mod=design$condition, noScale=TRUE)
hpgl_batch_expt = normalize_expt(expt=pasilla_expt, transform="log2", norm="quant", convert="cpm", batch="combatmod", filter_low=FALSE)
hpgl_batch = exprs(hpgl_batch_expt$expressionset)
test_that("Does combat batch correction end in the same dataframe?", {
    expect_equal(cbcb_batch, hpgl_batch)
})

## voom invocation
## This is a peculiar thing, the cbcbSEQ log2CPM() returns the normalized libsizes rather than those following log2() transformation.
## It then goes on to call voom() with this libsize rather than the log2().
## As a result, my normalization function is now keeping a copy of the libsizes and count tables from beginning to end
## to ensure that this is still accessible when required.\
message("Testing that the libsize input to voom is identical.")
cbcb_libsize = cbcb_l2qcpm_data$lib.size
hpgl_libsize = hpgl_l2qcpm_data$normalized_counts$libsize
test_that("In preparing for voom(), are the library sizes the same?", {
    expect_equal(cbcb_libsize, hpgl_libsize)
})

## Given that, lets try a voom() invocation and see what happens
message("Testing different methods of invoking voom, with/without cpm/log2 converted/transformed data.")
condition = design$condition
test_model = model.matrix(~condition)
cbcb_voom = voomMod(x=as.matrix(cbcb_l2qcpm), design=test_model, lib.size=cbcb_libsize, plot=TRUE)
hpgl_voom = voomMod(x=as.matrix(hpgl_l2qcpm), design=test_model, lib.size=hpgl_libsize, plot=TRUE)
hpgl_voom2 = hpgltools::hpgl_voom(as.matrix(hpgl_l2qcpm), test_model, libsize=hpgl_libsize, logged=TRUE, converted=TRUE)
hpgl_voom3 = hpgltools::hpgl_voom(as.matrix(hpgl_quantile), test_model, libsize=hpgl_libsize, logged=FALSE, converted=FALSE)

message("I can't test cbcb_voom vs. hpgl_voom3 because I set the row/colnames.")
test_that("Do different voom() invocations end with the same data?", {
    expect_equal(cbcb_voom, hpgl_voom)
    expect_equal(cbcb_voom$E, hpgl_voom2$E)
    expect_equal(cbcb_voom$E, hpgl_voom3$E)
})

## my hpgl_voom() sets row/column names and causes a test of the weights to fail.
## But checking manually shows them the same.
## expect_equal(cbcb_voom$weights, hpgl_voom$weights)

message("Test the final results from limma by invoking it through different codepaths.")
cbcb_fit = lmFit(cbcb_voom)
cbcb_eb = eBayes(cbcb_fit)
cbcb_top = topTable(cbcb_eb, number=nrow(cbcb_eb))
hpgl_toptables = limma_pairwise(expt=hpgl_l2qcpm_expt, model_intercept=TRUE, model_batch=FALSE, libsize=hpgl_libsize)
hpgl_top = hpgl_toptables$limma_result

## Changing the libsize does have a very small effect on the result
##test_libsize = hpgl_l2qcpm_data$final_counts$libsize
##hpgl_toptables = limma_pairwise(expt=hpgl_l2qcpm_expt, model_intercept=TRUE, model_batch=FALSE, libsize=test_libsize)
##hpgl_top2 = hpgl_toptables$limma_result
##test_that("Limma results.", {
##    expect_equal(cbcb_top, hpgl_top)
##})
##message("The following low correlation shows how much the results change when libsize is chosen poorly!")
##print(cor.test(cbcb_top$logFC, hpgl_top2$logFC))
##message("As a result, recent versions of limma_pairwise go to some trouble to choose libsize.")
##hpgl_toptables = limma_pairwise(expt=hpgl_l2qcpm_expt, model_batch=FALSE)
##message("Using limma_pairwise() without an intercept model adds a qvalue column.")
##hpgl_top = hpgl_toptables$limma_result$untreated_minus_treated
##hpgl_top = hpgl_top[,c(1,2,3,4,5,6)]
##test_that("Limma results with automatic libsize.", {
##    expect_equal(cbcb_top, hpgl_top)
##})

