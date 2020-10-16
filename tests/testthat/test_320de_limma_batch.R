start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
cbcb <- sm(library(cbcbSEQ))
context("320de_limma_batch.R: Does hpgltools work with limma?\n")

load("pasilla_df.rda")
pasilla <- new.env()
load("pasilla.rda", envir=pasilla)
pasilla_expt <- pasilla[["expt"]]

## Testing that hpgltools gets a similar result to cbcbSEQ using limma.
## Until preprocessCore gets fixed, disable this.
cbcb_counts <- cbcbSEQ::filterCounts(counts)
cbcb_qcounts <- cbcbSEQ::qNorm(cbcb_counts)
cbcb_cpm <- cbcbSEQ::log2CPM(cbcb_qcounts)
cbcb_qcpmcounts <- as.matrix(cbcb_cpm[["y"]])
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_qcpmcounts)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd[["v"]], cbcb_svd[["d"]],
                           design[["condition"]], design[["libType"]])

cbcb_libsize <- cbcb_cpm[["lib.size"]]
cbcb_combat <- cbcbSEQ::combatMod(cbcb_cpm[["y"]], batch=design[["libType"]],
                                  mod=design[["condition"]], noScale=TRUE)
## oh yeah, cbcbSEQ's combatMod no longer works
cbcb_v <- cbcbSEQ::voomMod(cbcb_qcpmcounts,
                           model.matrix(~ design[["condition"]] + design[["libType"]]),
                           lib.size=cbcb_libsize)
## It looks to me like the voomMod function is missing a is.na() check and so
## the lowess() function is failing.
hpgl_v <- hpgl_voom(cbcb_qcpmcounts,
                    model=model.matrix(~ design[["condition"]] + design[["libType"]]),
                    libsize=cbcb_libsize, logged=TRUE, converted=TRUE)
## Taking the first column of the E slot in in v

cbcb_fit <- lmFit(cbcb_v)
cbcb_eb <- eBayes(cbcb_fit)
table_order <- sort(rownames(cbcb_eb))
cbcb_table <- topTable(cbcb_eb, coef=2, n=nrow(cbcb_v[["E"]]))

cbcb_data <- as.matrix(counts)
cbcb_data <- cbcb_data[table_order, ]
hpgl_data <- exprs(pasilla_expt)
hpgl_data <- hpgl_data[table_order, ]
test_that("Does data from an expt equal a raw dataframe?", {
  expect_equal(cbcb_data, hpgl_data)
})

## Perform log2/cpm/quantile/combatMod normalization
hpgl_norm <- sm(normalize_expt(pasilla_expt, transform="log2", norm="quant",
                               convert="cbcbcpm", filter="cbcb", thresh=1))

## Ensure that we have the same count tables for limma_pairwise
## and the invocations of voom->topTable() by cbcbSEQ.
##expected <- nrow(cbcb_counts)
##actual <- nrow(exprs(hpgl_norm))
##test_that("Do we get the same number of genes using cbcb's filter as normalize_expt?", {
##  expect_equal(expected, actual)
##})

## If we made it this far, then the inputs to limma should agree.
## Use this section to ensure that our invocation of limma without an intercept
## matches that from cbcbSEQ.

int_limma <- sm(limma_pairwise(hpgl_norm, model_batch=TRUE, limma_method="ls",
                               model_intercept=TRUE, which_voom="hpgl"))
int_voom <- int_limma[["voom_result"]]
int_fit <- int_limma[["fit"]]
int_eb <- int_limma[["pairwise_comparisons"]]
int_table <- int_limma[["all_tables"]][["untreated"]]

## Now see that the voom outputs are the same
expected <- cbcb_v[["E"]]
expected <- expected[table_order, ]
actual <- int_voom[["E"]]
actual <- actual[table_order, ]
test_that("Do cbcbSEQ and hpgltools agree on the voom output?", {
  expect_equal(expected, actual)
})

## Fix the column names for the following tables to simplify things.
int_names <- c("(Intercept)", "untreated", "single_end")

expected <- cbcb_fit[["coefficients"]]
colnames(expected) <- int_names
expected <- expected[table_order, ]
actual <- int_fit[["coefficients"]]
actual <- actual[table_order, ]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: coefficients?", {
  expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["stdev.unscaled"]]
colnames(expected) <- int_names
expected <- expected[table_order, ]
actual <- int_fit[["stdev.unscaled"]]
actual <- actual[table_order, ]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: stdev_unscaled?", {
  expect_equal(expected, actual, tolerance=0.001)
})

expected <- cbcb_fit[["cov.coefficients"]]
colnames(expected) <- int_names
rownames(expected) <- int_names
actual <- int_fit[["cov.coefficients"]]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result: cov.coefficients?", {
  expect_equal(expected, actual)
})

expected <- cbcb_eb[["t"]]
colnames(expected) <- int_names
expected <- expected[table_order, ]
actual <- int_eb[["t"]]
actual <- actual[table_order, ]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result: eb[1]?", {
  expect_equal(expected, actual, tolerance=0.1) ## The intercept
})

expected <- cbcb_eb[["p.value"]]
colnames(expected) <- int_names
expected <- expected[table_order, ]
actual <- int_eb[["p.value"]]
actual <- actual[table_order, ]
test_that("Do the p-value tables stay the same pval[1]?", {
  expect_equal(expected[[1]], actual[[1]])
})

expected <- cbcb_table[table_order, "logFC"]
actual <- int_table[table_order, "logFC"]
test_that("Do cbcbSEQ and hpgltools agree on the logFCs?", {
  expect_equal(expected, actual, tolerance=0.01)
})

expected <- cbcb_table[table_order, "AveExpr"]
actual <- int_table[table_order, "AveExpr"]
test_that("Do cbcbSEQ and hpgltools agree on the AveExprs?", {
  expect_equal(expected, actual, tolerance=0.01)
})

expected <- cbcb_table[table_order, "P.Value"]
actual <- int_table[table_order, "P.Value"]
test_that("Do cbcbSEQ and hpgltools agree on the p-values?", {
  expect_equal(expected, actual, tolerance=0.01)
})

expected <- cbcb_table[table_order, "adj.P.Value"]
actual <- int_table[table_order, "adj.P.Value"]
test_that("Do cbcbSEQ and hpgltools agree on the p-values?", {
  expect_equal(expected, actual, tolerance=0.01)
})

## Finished checking the no-intercept invocations, now compare intercept to no-intercept limma.
noint_limma <- sm(limma_pairwise(hpgl_norm, which_voom="hpgl", limma_method="ls"))
expected <- noint_limma[["voom_result"]][["E"]]
table_order <- rownames(expected)
actual <- int_limma[["voom_result"]][["E"]]
expected <- expected[table_order, ]
actual <- actual[table_order, ]
test_that("Are the intercept and non-intercept voom results equivalent?", {
  expect_equal(expected, actual)
})

noint_coefficients <- noint_limma[["fit"]][["coefficients"]]
int_coefficients <- int_limma[["fit"]][["coefficients"]]
expected <- int_fit[["coefficients"]][, 1]  ## The '(Intercept)' column
actual <- noint_limma[["fit"]][["coefficients"]][, "treated"]
test_that("Do the intercept model results equal those from cell means for the intercept?", {
  expect_equal(expected, actual)
})

## Now check the noint_table vs. the int_table results.
noint_table <- noint_limma[["all_tables"]][[1]]
noint_table <- noint_table[table_order, ]
int_table <- int_limma[["all_tables"]][[1]]
int_table <- int_table[table_order, ]
expected <- noint_table[, "logFC"]
actual <- int_table[, "logFC"]
test_that("Do the intercept model results equal those from no-intercept (logFC)?", {
  expect_equal(expected, actual)
})

expected <- noint_table[["AveExpr"]]
actual <- int_table[["AveExpr"]]
test_that("Do the intercept and no-intercept fits give equal AveExpr values?", {
  expect_equal(expected, actual)
})

expected <- noint_table[["t"]]
actual <- int_table[["t"]]
test_that("Do the intercept and no-intercept fits give equal t-statistics?", {
  expect_equal(expected, actual, tolerance=0.02)
})

expected <- noint_table[["P.Value"]]
actual <- int_table[["P.Value"]]
test_that("Do the intercept and no-intercept fits give equal P-Values?", {
  expect_equal(expected, actual, tolerance=0.01)
})

limma_written <- sm(write_limma(noint_limma, excel="limma_test.xlsx"))
hpgl_limma <- sm(limma_pairwise(pasilla_expt))

## For the following tests
limma_file <- "320_de_limma.rda"
saved <- save(list=ls(), file=limma_file)
test_that("Did we save the limma results?", {
  expect_true(file.exists(limma_file))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 320de_limma_batch.R in ", elapsed,  " seconds."))
