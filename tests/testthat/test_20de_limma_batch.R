start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
cbcb <- sm(library(cbcbSEQ))
context("20de_limma_batch.R: Does limma work with hpgltools?\n")

load("pasilla_df.rda")

pasilla <- new.env()
load("pasilla.Rdata", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## Testing that hpgltools gets a similar result to cbcbSEQ using limma.
cbcb_counts <- cbcbSEQ::filterCounts(counts)
cbcb_qcounts <- cbcbSEQ::qNorm(cbcb_counts)
cbcb_cpm <- cbcbSEQ::log2CPM(cbcb_qcounts)
cbcb_qcpmcounts <- as.matrix(cbcb_cpm[["y"]])
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_qcpmcounts)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd[["v"]], cbcb_svd[["d"]],
                           design[["condition"]], design[["libType"]])

cbcb_libsize <- cbcb_cpm[["lib.size"]]
## cbcb_combat <- cbcbSEQ::combatMod(cbcb_cpm, batch=design[["libType"]],
##                                   mod=design[["condition"]], noScale=TRUE)
## oh yeah, cbcbSEQ's combatMod no longer works
cbcb_v <- cbcbSEQ::voomMod(cbcb_qcpmcounts,
                           model.matrix(~design[["condition"]] + design[["libType"]]),
                           lib.size=cbcb_libsize)
## It looks to me like the voomMod function is missing a is.na() check and so
## the lowess() function is failing.
hpgl_v <- hpgl_voom(cbcb_qcpmcounts,
                    model=model.matrix(~design[["condition"]] + design[["libType"]]),
                    libsize=cbcb_libsize, logged=TRUE, converted=TRUE)
## Taking the first column of the E slot in in v
cbcb_fit <- lmFit(cbcb_v)
cbcb_eb <- eBayes(cbcb_fit)
cbcb_table <- topTable(cbcb_eb, coef=2, n=nrow(cbcb_v[["E"]]))

cbcb_data <- as.matrix(counts)
cbcb_data <- cbcb_data[sort(rownames(cbcb_data)), ]
hpgl_data <- Biobase::exprs(pasilla_expt[["expressionset"]])
hpgl_data <- hpgl_data[sort(rownames(hpgl_data)), ]
test_that("Does data from an expt equal a raw dataframe?", {
    expect_equal(cbcb_data, hpgl_data)
})

## Perform log2/cpm/quantile/combatMod normalization
hpgl_norm <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant",
                               convert="cbcbcpm", filter=TRUE))

## If we made it this far, then the inputs to limma should agree.
hpgl_limma_nointercept <- sm(limma_pairwise(hpgl_norm, model_batch=TRUE, limma_method="ls",
                                            model_intercept=FALSE, which_voom="hpgl"))
hpgl_voom <- hpgl_limma_nointercept[["voom_result"]]
hpgl_fit <- hpgl_limma_nointercept[["fit"]]
hpgl_eb <- hpgl_limma_nointercept[["pairwise_comparisons"]]
hpgl_table <- hpgl_limma_nointercept[["all_tables"]][[1]]

hpgl_limma <- sm(limma_pairwise(hpgl_norm, which_voom="hpgl", limma_method="ls"))

expected <- cbcb_v[["E"]]
expected <- expected[sort(rownames(expected)), ]
actual <- hpgl_v[["E"]]
actual <- actual[sort(rownames(actual)), ]
test_that("Do cbcbSEQ and hpgltools agree on the voom output?", {
    expect_equal(expected, actual)
})

expected <- cbcb_fit[["coefficients"]]
expected <- as.numeric(head(expected[sort(rownames(expected)), ][[1]]))
actual <- hpgl_fit[["coefficients"]]
actual <- as.numeric(head(actual[sort(rownames(actual)), ][[1]]))
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: coefficients[1]?", {
    expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["coefficients"]]
expected <- as.numeric(head(expected[sort(rownames(expected)), ][[2]]))
actual <- hpgl_fit[["coefficients"]]
actual <- as.numeric(head(actual[sort(rownames(actual)), ][[2]]))
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: coefficients[2]?", {
    expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["stdev.unscaled"]]
expected <- as.numeric(head(expected[sort(rownames(expected)), ][[1]]))
actual <- hpgl_fit[["stdev.unscaled"]]
actual <- as.numeric(head(actual[sort(rownames(actual)), ][[1]]))
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: stdev_unscaled[1]?", {
    expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["stdev.unscaled"]]
expected <- as.numeric(head(expected[sort(rownames(expected)), ][[2]]))
actual <- hpgl_fit[["stdev.unscaled"]]
actual <- as.numeric(head(actual[sort(rownames(actual)), ][[2]]))
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: stdev_unscaled[2]?", {
    expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["cov.coefficients"]]
expected <- as.numeric(head(expected[sort(rownames(expected)), ][[1]]))
actual <- hpgl_fit[["cov.coefficients"]]
actual <- as.numeric(head(actual[sort(rownames(actual)), ][[1]]))
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: cov.coefficients?", {
    expect_equal(expected, actual)
})

expected <- cbcb_eb[["t"]]
expected <- expected[sort(rownames(expected)), ]
actual <- hpgl_eb[["t"]]
actual <- actual[sort(rownames(actual)), ]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result: eb[1]?", {
    expect_equal(expected[[1]], actual[[1]], tolerance=0.1) ## The intercept
})
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result: eb[2]?", {
    expect_equal(expected[[2]], actual[[2]], tolerance=0.1) ## condition-untreated
})
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result: eb[3]?", {
    expect_equal(expected[[3]], actual[[3]], tolerance=1) ## batch-single_end
})

expected <- cbcb_eb[["p.value"]]
expected <- expected[sort(rownames(expected)), ]
actual <- hpgl_eb[["p.value"]]
actual <- actual[sort(rownames(actual)), ]
test_that("Do the p-value tables stay the same pval[1]?", {
    expect_equal(expected[[1]], actual[[1]])
})
test_that("Do the p-value tables stay the same pval[2]?", {
    expect_equal(expected[[2]], actual[[2]], tolerance=0.01)
})
test_that("Do the p-value tables stay the same pval[3]?", {
    expect_equal(expected[[3]], actual[[3]])
})

## I need to be smarter about bringing together the filtering, otherwise these tests
## will continue to fail in weird ways.
cbcb_result_reordered <- cbcb_table[sort(rownames(actual)), ]
hpgl_result_reordered <- hpgl_table[sort(rownames(actual)), ]
cbcb_logfc <- as.numeric(head(cbcb_result_reordered[["logFC"]]))
hpgl_logfc <- as.numeric(head(hpgl_result_reordered[["untreated"]]))
test_that("Do cbcbSEQ and hpgltools agree on the list of DE genes?", {
    expect_equal(cbcb_logfc, hpgl_logfc, tolerance=0.0001)
})

reordered <- hpgl_limma[["all_tables"]][["untreated_vs_treated"]]
reordered <- reordered[sort(rownames(actual)), ]
test_that("Do the intercept model results equal those from cell means?", {
    expect_equal(hpgl_voom[["E"]], hpgl_limma[["voom_result"]][["E"]])
})
test_that("Do the intercept model results equal those from cell means?", {
    expect_equal(hpgl_fit[["coefficients"]][[1]], hpgl_limma[["fit"]][["coefficients"]][[1]])
})
## Something is not right here.
##test_that("Do the intercept model results equal those from cell means?", {
##    expect_equal(head(hpgl_eb[["p.value"]][[2]]),
##                 head(hpgl_limma[["pairwise_comparisons"]][["p.value"]][[1]]),
##                 tolerance=0.1)
##})
test_that("Do the intercept model results equal those from cell means?", {
    expect_equal(as.numeric(head(hpgl_logfc)), as.numeric(head(reordered[["logFC"]])), tolerance=0.1)
})

limma_written <- sm(write_limma(hpgl_limma, excel="limma_test.xlsx"))


save(list=ls(), file="de_limma.rda")

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 20de_limma_batch.R in ", elapsed,  " seconds."))
tt <- clear_session()
