library(testthat)
library(hpgltools)
context("Test usability of limma")
require.auto("pasilla")
require.auto("edgeR")
require.auto("kokrah/cbcbSEQ")

context("Testing limma usage.")

message("Loading pasilla, setting up count tables.")
message("Taking this section directly from the cbcbSEQ vignette.")
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
counts = counts[rowSums(counts) > ncol(counts), ]
design = data.frame(row.names=colnames(counts),
    condition=c("untreated","untreated","untreated",
        "untreated","treated","treated","treated"),
    libType=c("single-end","single-end","paired-end",
        "paired-end","single-end","paired-end","paired-end"))
metadata = design
colnames(metadata) = c("condition", "batch")
metadata$Sample.id = rownames(metadata)

message("Testing that hpgltools gets a similar result to cbcbSEQ using limma.")
library(cbcbSEQ)
cbcb_qcounts <- cbcbSEQ::qNorm(counts)
cbcb_cpm <- cbcbSEQ::log2CPM(cbcb_qcounts)
cbcb_qcpmcounts <- as.matrix(cbcb_cpm[["y"]])
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_qcpmcounts)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd$v, cbcb_svd$d, design$condition, design$libType)
cbcb_vignette_result <- c(27.57, 24.66, 15.62, 12.15, 10.53, 9.46)
test_that("Does cbcbSEQ give the same result for the initial pcRes call?", {
    expect_equal(cbcb_vignette_result, as.numeric(cbcb_res$propVar))
})

cbcb_libsize <- cbcb_cpm[["lib.size"]]
## cbcb_combat <- cbcbSEQ::combatMod(cbcb_cpm, batch=design[["libType"]], mod=design[["condition"]], noScale=TRUE)
## oh yeah, cbcbSEQ's combatMod no longer works
cbcb_hpgl_combat <- hpgl_combatMod(dat=cbcb_qcpmcounts, batch=design[["libType"]], mod=design[["condition"]], noScale=TRUE)
## Ok, here is a point where the cbcbSEQ vignette does not agree with its output.
## the return of cbcbSEQ::combatMod (if it worked) is a variable containing only 'bayesdata', not a list of bayesdata and info.
message("Test again that cbcbSEQ's principle components match these (since I dropped in a different implementation of combatMod.")
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_hpgl_combat)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd$v, cbcb_svd$d, design$condition, design$libType)
cbcb_almost_vignette_result <- c(30.39, 18.56, 14.71, 12.92, 12.39, 11.03)  ## Taken from when I run the commands in the vignette.
cbcb_actual_vignette_result <- c(30.97, 18.65, 14.69, 12.65, 12.09, 10.94)  ## Taken from cbcbSEQIntro.pdf
test_that("Does the post-batch correction PCA give the same result?", {
    expect_equal(cbcb_almost_vignette_result, as.numeric(cbcb_res$propVar))
})
cbcb_v <- cbcbSEQ::voomMod(cbcb_hpgl_combat, model.matrix(~design$condition), lib.size=cbcb_libsize, plot=FALSE)
## It looks to me like the voomMod function is missing a is.na() check and so the lowess() function is failing.
hpgl_v <- hpgl_voom(cbcb_hpgl_combat, model=model.matrix(~design$condition), libsize=cbcb_libsize, logged=TRUE, converted=TRUE)
## Taking the first column of the E slot in in v
cbcb_almost_vignette_result <- c(2.968411, 3.028748, 3.265501, 2.858357, 2.838402, 3.178890, 2.713208)
cbcb_actual_vignette_result <- c(2.9772407, 3.0375781, 3.259578, 2.852434, 2.847232, 3.1729673, 2.7072849)
test_that("Does the cbcbSEQ voomMod() function give the same results as hpgl_voom()?", {
    expect_equal(cbcb_v$E, hpgl_v$E) })
test_that("Do they agree with my approximated vignette results?", {
    expect_equal(as.numeric(head(cbcb_v$E, n=1)), cbcb_almost_vignette_result, tolerance=0.0001)
})
cbcb_fit <- lmFit(cbcb_v)
cbcb_eb <- eBayes(cbcb_fit)
cbcb_table <- topTable(cbcb_eb, coef=2, n=nrow(v$E))


## Now create a hpgltools expt and try the same thing
message("Setting up an expt class to contain the pasilla data and metadata.")
pasilla_expt = create_expt(count_dataframe=counts, meta_dataframe=metadata)
cbcb_data = as.matrix(counts)
hpgl_data = Biobase::exprs(pasilla_expt$expressionset)
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(cbcb_data, hpgl_data)
})

## Perform log2/cpm/quantile/combatMod normalization
hpgl_norm <- normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm")
hpgl_qcpmcounts <- Biobase::exprs(hpgl_norm$expressionset)
hpgl_qcpm_combat_counts <- normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm", batch="combatmod")
test_that("Do cbcbSEQ and hpgltools agree on the definition of log2(quantile(cpm(counts)))?", {
    expect_equal(cbcb_qcpmcounts, hpgl_qcpmcounts)
})

hpgl_qcpmcombat <- normalize_expt(pasilla_expt, transform="log2", norm="quant", convert="cpm", batch="combatmod")
hpgl_combat <- Biobase::exprs(hpgl_qcpmcombat$expressionset)
test_that("Do cbcbSEQ and hpgltools agree on combatMod(log2(quantile(cpm(counts))))?", {
    expect_equal(cbcb_hpgl_combat, hpgl_combat)
})

message("If we made it this far, then the inputs to limma should agree.")
hpgl_limma_result <- limma_pairwise(hpgl_qcpmcombat, model_batch=FALSE, model_intercept=TRUE)
hpgl_voom <- hpgl_limma_result$voom_result
hpgl_fit <- hpgl_limma_result$fit
hpgl_eb <- hpgl_limma_result$pairwise_comparisons
hpgl_table <- hpgl_limma_result$all_tables
message("The order of operations in a limma analysis are: voom->fit->ebayes->table, test them in that order.")
message("Keep in mind that I do not default to an intercept model, and I rename the columns of the coefficients to make them more readable.")
test_that("Do cbcbSEQ and hpgltools agree on the voom output?", {
    expect_equal(cbcb_v$E, hpgl_voom$E)
})
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result?", {
    expect_equal(cbcb_fit$coefficients[[1]], hpgl_fit$coefficients[[1]])
    expect_equal(cbcb_fit$coefficients[[2]], hpgl_fit$coefficients[[2]])
    expect_equal(cbcb_fit$stdev.unscaled[[1]], hpgl_fit$stdev.unscaled[[1]])
    expect_equal(cbcb_fit$stdev.unscaled[[2]], hpgl_fit$stdev.unscaled[[2]])
    expect_equal(cbcb_fit$df.residual, hpgl_fit$df.residual)
    expect_equal(cbcb_fit$cov.coefficients[[1]], hpgl_fit$cov.coefficients[[1]])
    expect_equal(cbcb_fit$cov.coefficients[[2]], hpgl_fit$cov.coefficients[[2]])
    expect_equal(cbcb_fit$pivot, hpgl_fit$pivot)
    expect_equal(cbcb_fit$rank, hpgl_fit$rank)
    expect_equal(cbcb_fit$Amean, hpgl_fit$Amean)
    message("There are a couple more slots in these returns, but I think by now we get the point.")
})

test_that("Do cbcbSEQ and hpgltools agree on the eBayes result?", {
    expect_equal(cbcb_eb$t[[1]], hpgl_eb$t[[1]])
    expect_equal(cbcb_eb$t[[2]], hpgl_eb$t[[2]])
    expect_equal(cbcb_eb$p.value[[1]], hpgl_eb$p.value[[1]])
    expect_equal(cbcb_eb$p.value[[2]], hpgl_eb$p.value[[2]])
    message("The eBayes results include the previous fits and some more slots.  I only tested a few here.")
})

test_that("Do cbcbSEQ and hpgltools agree on the list of DE genes?", {
    expect_equal(cbcb_table, hpgl_table)
})

message("YAY!")
