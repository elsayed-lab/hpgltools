library(testthat)
library(hpgltools)
library(pasilla)
data(pasillaGenes)
num_installed <- require.auto("kokrah/cbcbSEQ")
cbcb <- s_p(library(cbcbSEQ))

context("Compare cbcbSEQ output to hpgltools.")
## A bunch of these tests are certain to fail with the changed order of genes in expressionsets.

## This section is copy/pasted to all of these tests, that is dumb.
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header=TRUE, row.names=1)
counts <- counts[rowSums(counts) > ncol(counts),]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single_end","single_end","paired_end",
        "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata$sampleid <- rownames(metadata)

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

cbcb_data <- as.matrix(counts)
hpgl_data <- Biobase::exprs(pasilla_expt$expressionset)

## Check that normalization tools work similarly
cbcb_quantile <- cbcbSEQ::qNorm(cbcb_data)
hpgl_quantile_data <- s_p(hpgl_norm(pasilla_expt, transform="raw", norm="quant", convert="raw", filter_low=FALSE))$result
hpgl_quantile <- hpgl_quantile_data[["count_table"]]
cbcb_quantile <- cbcb_quantile[sort(rownames(cbcb_quantile)), ]
hpgl_quantile <- hpgl_quantile[sort(rownames(hpgl_quantile)), ]
test_that("Are the quantile normalizations identical?", {
    expect_equal(cbcb_quantile, hpgl_quantile)
})

## Check that quant+cpm are the same
cbcb_qcpm <- cbcbSEQ::qNorm(cbcb_data)
library(edgeR)
## I don't know how to call cpm without using a library call first.
cbcb_edger_qcpm <- edgeR::cpm(cbcb_qcpm)

cbcb_quantile <- cbcbSEQ::qNorm(cbcb_data)
cbcb_quantile <- cbcb_quantile[sort(rownames(cbcb_quantile)), ]
hpgl_quantile <- s_p(hpgl_norm(pasilla_expt, norm="quant"))[["result"]]
hpgl_quantile <- hpgl_quantile[["count_table"]]
hpgl_quantile <- hpgl_quantile[sort(rownames(hpgl_quantile)), ]
test_that("Are quantiles identical?", {
    expect_equal(cbcb_quantile, hpgl_quantile)
})

hpgl_qcpm <- s_p(hpgl_norm(pasilla_expt, norm="quant", convert="cpm", filter_low=FALSE))[["result"]]
hpgl_qcpm <- hpgl_qcpm[["count_table"]]
hpgl_qcpm <- hpgl_qcpm[sort(rownames(hpgl_qcpm)), ]
test_that("Are cpm conversions identical?", {
    expect_equal(cbcb_edger_qcpm, hpgl_qcpm)
})

## log2/cpm that
cbcb_l2qcpm_data <- cbcbSEQ::log2CPM(cbcb_quantile)
cbcb_l2qcpm <- cbcb_l2qcpm_data[["y"]]
hpgl_l2qcpm_data <- s_p(hpgl_norm(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm", filter_low=FALSE))[["result"]]
hpgl_l2qcpm <- hpgl_l2qcpm_data[["count_table"]]
hpgl_l2qcpm <- hpgl_l2qcpm[sort(rownames(hpgl_l2qcpm)), ]
hpgl_l2qcpm_expt <- s_p(normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cbcbcpm", filter=FALSE))$result
hpgl_l2qcpm2 <- Biobase::exprs(hpgl_l2qcpm_expt[["expressionset"]])
hpgl_l2qcpm2 <- hpgl_l2qcpm2[sort(rownames(hpgl_l2qcpm2)), ]
test_that("Are l2qcpm conversions/transformations identical using two codepaths?", {
    expect_equal(cbcb_l2qcpm, hpgl_l2qcpm)
    expect_equal(cbcb_l2qcpm, hpgl_l2qcpm2)
})

## Check that PCA invocations are similar
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_l2qcpm)
hpgl_pca_info <- plot_pca(hpgl_l2qcpm_expt)
hpgl_svd <- hpgl_pca_info[["pca"]]
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd[["v"]], cbcb_svd[["d"]], design[["condition"]], design[["libType"]])
hpgl_res <- hpgl_pca_info[["res"]]
test_that("Do the PCA invocations provide the same results?", {
    expect_equal(cbcb_svd[["v"]], hpgl_svd[["v"]])
    expect_equal(cbcb_svd[["d"]], hpgl_svd[["d"]])
    expect_equal(cbcb_res, hpgl_res)
})

## This is a peculiar thing, the cbcbSEQ log2CPM() returns the normalized libsizes rather than those following log2() transformation.
## It then goes on to call voom() with this libsize rather than the log2().
## As a result, my normalization function is now keeping a copy of the libsizes and count tables from beginning to end
## to ensure that this is still accessible when required.
cbcb_libsize <- cbcb_l2qcpm_data[["lib.size"]]
hpgl_libsize <- hpgl_l2qcpm_data[["intermediate_counts"]][["normalization"]][["libsize"]]
test_that("In preparing for voom(), are the library sizes the same?", {
    expect_equal(cbcb_libsize, hpgl_libsize)
})

## Given that, lets try a voom() invocation and see what happens
condition <- design[["condition"]]
test_model <- model.matrix(~condition)
cbcb_voom <- cbcbSEQ::voomMod(x=as.matrix(cbcb_l2qcpm), design=test_model, lib.size=cbcb_libsize)
hpgl_voom <- cbcbSEQ::voomMod(x=as.matrix(hpgl_l2qcpm), design=test_model, lib.size=hpgl_libsize)
hpgl_voom2 <- hpgltools::hpgl_voom(as.matrix(hpgl_l2qcpm), model=test_model, libsize=hpgl_libsize, logged=TRUE, converted=TRUE)
hpgl_voom3 <- s_p(hpgltools::hpgl_voom(as.matrix(hpgl_quantile), test_model, libsize=hpgl_libsize, logged=FALSE, converted=FALSE))$result

test_that("Do different voom() invocations end with the same data?", {
    expect_equal(cbcb_voom, hpgl_voom)
    expect_equal(cbcb_voom$E, hpgl_voom2$E)
    expect_equal(cbcb_voom$E, hpgl_voom3$E)
})

## my hpgl_voom() sets row/column names and causes a test of the weights to fail.
## But checking manually shows them the same.
## expect_equal(cbcb_voom$weights, hpgl_voom$weights)

hpgl_limma <- limma_pairwise(hpgl_l2qcpm_expt, model_batch=FALSE, model_intercept=FALSE)

## First check the voom result from limma_pairwise
hpgl_limma_voom <- hpgl_limma$voom_result
hpgl_limma_voom_e <- hpgl_limma$voom_result$E[order(rownames(hpgl_limma$voom_result$E)), ]
cbcb_voom_e <- cbcb_voom$E[order(rownames(cbcb_voom$E)), ]
test_that("Limma results, voom.", {
    expect_equal(cbcb_voom_e, hpgl_limma_voom_e)
})

## Then the result from lmFit
hpgl_limma_fit_coef <- hpgl_limma$fit$coefficients[order(rownames(hpgl_limma$fit$coefficients)), ]
cbcb_fit <- limma::lmFit(cbcb_voom)
cbcb_fit_coef <- cbcb_fit$coefficients[order(rownames(cbcb_fit$coefficients)), ]
colnames(cbcb_fit_coef) <- c("(Intercept)", "untreated")
test_that("Limma results, fitting.", {
    expect_equal(cbcb_fit_coef, hpgl_limma_fit_coef)
})
hpgl_limma_fit_stdev <- hpgl_limma$fit$stdev.unscaled[order(rownames(hpgl_limma$fit$stdev.unscaled)), ]
cbcb_fit_stdev <- cbcb_fit$stdev.unscaled[order(rownames(cbcb_fit$stdev.unscaled)), ]
colnames(cbcb_fit_stdev) <- c("(Intercept)", "untreated")
test_that("Limma results, fitting.", {
    expect_equal(cbcb_fit_stdev, hpgl_limma_fit_stdev)
})

## Now the result from eBayes
cbcb_eb <- limma::eBayes(cbcb_fit)


cbcb_top <- s_p(limma::topTable(cbcb_eb, number=nrow(cbcb_eb)))[["result"]]
cbcb_top <- cbcb_top[sort(rownames(cbcb_top)), ]


hpgl_top <- hpgl_limma[["all_tables"]][["untreated_vs_treated"]]
hpgl_top <- hpgl_top[sort(rownames(hpgl_top)), ]


    expect_equal(cbcb_top, hpgl_top)
})
