start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("321de_limma_combat.R: Does hpgltools work with limma and combat?")

## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]
limma <- new.env()
load("320_de_limma.rda", envir = limma)
counts <- limma[["counts"]]
design <- limma[["design"]]

## Testing that hpgltools gets a similar result to cbcbSEQ using limma.
cbcb <- sm(library(cbcbSEQ))
counts <- limma[["counts"]]
design <- limma[["design"]]
cbcb_qcounts <- cbcbSEQ::qNorm(counts)
cbcb_cpm <- cbcbSEQ::log2CPM(cbcb_qcounts)
cbcb_qcpmcounts <- as.matrix(cbcb_cpm[["y"]])
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_qcpmcounts)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd[["v"]], cbcb_svd[["d"]],
                           design[["condition"]], design[["libType"]])
cbcb_vignette_result <- c(27.57, 24.66, 15.62, 12.15, 10.53, 9.46)
test_that("Does cbcbSEQ give the same result for the initial pcRes call?", {
  expect_equal(cbcb_vignette_result, as.numeric(cbcb_res[["propVar"]]))
})

cbcb_libsize <- cbcb_cpm[["lib.size"]]
## cbcb_combat <- cbcbSEQ::combatMod(cbcb_cpm, batch = design[["libType"]],
##                                   mod = design[["condition"]], noScale = TRUE)
## oh yeah, cbcbSEQ's combatMod no longer works
cbcb_hpgl_combat <- sm(cbcb_combat(dat = cbcb_qcpmcounts,
                                   batch = design[["libType"]],
                                   mod = design[["condition"]], noScale = TRUE))
## Ok, here is a point where the cbcbSEQ vignette does not agree with its output.
## the return of cbcbSEQ::combatMod (if it worked) is a variable containing only
## 'bayesdata', not a list of bayesdata and info.

## Test again that cbcbSEQ's principle components match these (since I dropped
## in a different implementation of combatMod.
cbcb_svd <- cbcbSEQ::makeSVD(cbcb_hpgl_combat)
cbcb_res <- cbcbSEQ::pcRes(cbcb_svd[["v"]], cbcb_svd[["d"]],
                           design[["condition"]], design[["libType"]])

## The following was taken from when I run the commands in the vignette.
cbcb_almost_vignette_result <- c(30.39, 18.56, 14.71, 12.92, 12.39, 11.03)
## The following was taken from the cbcbSEQIntro.pdf
cbcb_actual_vignette_result <- c(30.97, 18.65, 14.69, 12.65, 12.09, 10.94)
test_that("Does the post-batch correction PCA give the same result?", {
  expect_equal(cbcb_almost_vignette_result, as.numeric(cbcb_res[["propVar"]]))
})
cbcb_v <- suppressWarnings(cbcbSEQ::voomMod(cbcb_hpgl_combat,
                                            model.matrix(~design[["condition"]]),
                                            lib.size = cbcb_libsize))
## It looks to me like the voomMod function is missing a is.na() check and so
## the lowess() function is failing.
hpgl_v <- suppressWarnings(hpgl_voom(cbcb_hpgl_combat,
                                     model = model.matrix(~design[["condition"]]),
                                     libsize = cbcb_libsize,
                                     logged = TRUE, converted = TRUE))
## Taking the first column of the E slot in in v
cbcb_almost_vignette_result <- c(2.968411, 3.028748, 3.265501, 2.858357,
                                 2.838402, 3.178890, 2.713208)
cbcb_actual_vignette_result <- c(2.9772407, 3.0375781, 3.259578, 2.852434,
                                 2.847232, 3.1729673, 2.7072849)
test_that("Does the cbcbSEQ voomMod() function give the same results as hpgl_voom()?", {
  expect_equal(cbcb_v[["E"]], hpgl_v[["E"]])
})
test_that("Do they agree with my approximated vignette results?", {
  expect_equal(as.numeric(head(cbcb_v[["E"]], n = 1)),
               cbcb_almost_vignette_result, tolerance = 0.0001)
})
cbcb_fit <- lmFit(cbcb_v)
cbcb_eb <- eBayes(cbcb_fit)
cbcb_table <- topTable(cbcb_eb, coef = 2, n = nrow(cbcb_v[["E"]]))

## Now create a hpgltools expt and try the same thing
## Setting up an expt class to contain the pasilla data and metadata.
expected <- as.matrix(counts)
expected <- expected[sort(rownames(expected)), ]
actual <- exprs(pasilla_expt)
actual <- actual[sort(rownames(actual)), ]
test_that("Does data from an expt equal a raw dataframe?", {
  expect_equal(expected, actual)
})

## Perform log2/cpm/quantile/combatMod normalization
hpgl_norm <- sm(normalize_expt(pasilla_expt, transform = "log2",
                               norm = "quant", convert = "cbcbcpm"))
hpgl_qcpmcounts <- exprs(hpgl_norm)
expected <- cbcb_qcpmcounts
expected <- expected[sort(rownames(expected)), ]
actual <- hpgl_qcpmcounts
actual <- actual[sort(rownames(actual)), ]
test_that("Do cbcbSEQ and hpgltools agree on the definition of log2(quantile(cpm(counts)))?", {
  expect_equal(expected, actual)
})

## Getting log2(combat(cpm(quantile(counts))))
cbcb_hpgl_combat <- cbcbSEQ::combatMod(dat = cbcb_qcpmcounts, batch = design[["libType"]],
                                       mod = design[["condition"]], noScale = TRUE)
cbcb_filtered <- suppressWarnings(cbcbSEQ::filterCounts(cbcb_hpgl_combat))
not_na <- !is.na(rownames(cbcb_filtered))
cbcb_filtered <- cbcb_filtered[not_na, ]

intermediate <- normalize_expt(pasilla_expt, filter = "cbcb", thresh = 2)
end <- normalize_expt(intermediate, norm = "quant", convert = "cbcbcpm", batch = "combatmod")
end <- normalize_expt(end, transform = "log2")

expected <- head(cbcb_filtered)[, 1]
actual <- head(exprs(end))[, 1]
test_that("Do cbcbSEQ and hpgltools agree on combatMod(log2(quantile(cpm(counts))))?", {
  expect_equal(expected, actual, tolerance = 0.1)
})

## If we made it this far, then the inputs to limma should agree.
hpgl_limma_combat_result <- sm(limma_pairwise(pasilla_expt, limma_method = "ls", filter = TRUE,
                                              model_batch = FALSE,
                                              model_intercept = TRUE,
                                              which_voom = "hpgl"))
hpgl_voom <- hpgl_limma_combat_result[["voom_result"]]
hpgl_fit <- hpgl_limma_combat_result[["fit"]]
hpgl_eb <- hpgl_limma_combat_result[["pairwise_comparisons"]]
hpgl_table <- hpgl_limma_combat_result[["all_tables"]][["untreated"]]

chosen_genes <- head(rownames(cbcb_filtered))
## The order of operations in a limma analysis are: voom->fit->ebayes->table,
## test them in that order.  Keep in mind that I do not default to an intercept
## model, and I rename the columns of the coefficients to make them more
## readable.
expected <- cbcb_v[["E"]][chosen_genes, 1]
actual <- hpgl_voom[["E"]][chosen_genes, 1]
test_that("Do cbcbSEQ and hpgltools agree on the voom output?", {
  expect_equal(expected, actual, tolerance = 0.1)
})

expected <- cbcb_fit[["coefficients"]][chosen_genes, 1]
actual <- hpgl_fit[["coefficients"]][chosen_genes, 1]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result? (Intercept)", {
  expect_equal(expected, actual, tolerance = 0.1)
})

expected <- cbcb_fit[["coefficients"]][chosen_genes, 2]
actual <- hpgl_fit[["coefficients"]][chosen_genes, 2]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result? (untreated)", {
  expect_equal(expected, actual, tolerance = 0.2)
})

expected <- cbcb_fit[["stdev.unscaled"]][chosen_genes, 1]
actual <- hpgl_fit[["stdev.unscaled"]][chosen_genes, 1]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result? (stdev.unscaled, intercept)", {
  expect_equal(expected, actual, tolerance = 0.1)
})

expected <- cbcb_fit[["stdev.unscaled"]][chosen_genes, 2]
actual <- hpgl_fit[["stdev.unscaled"]][chosen_genes, 2]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result? (stdev.unscaled, untreated)", {
  expect_equal(expected, actual, tolerance = 0.1)
})

expected <- cbcb_fit[["Amean"]][chosen_genes]
actual <- hpgl_fit[["Amean"]][chosen_genes]
test_that("Do cbcbSEQ and hpgltools agree on the lmFit result? (Amean)", {
  expect_equal(expected, actual, tolerance = 0.1)
})

expected <- cbcb_eb[["t"]]
expected <- expected[sort(rownames(expected)), ]
actual <- hpgl_eb[["t"]]
actual <- actual[sort(rownames(actual)), ]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result? (t, intercept)", {
  expect_equal(expected[[1]], actual[[1]], tolerance = 0.2)
})

expected <- cbcb_eb[["t"]][chosen_genes, 2]
actual <- cbcb_eb[["t"]][chosen_genes, 2]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result? (t, untreated)", {
  expect_equal(expected, actual, tolerance = 0.02)
})

expected <- cbcb_eb[["p.value"]][chosen_genes, 1]
actual <- hpgl_eb[["p.value"]][chosen_genes, 1]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result? (p.value, intercept)", {
  expect_equal(expected, actual, tolerance = 0.001)
})

expected <- cbcb_eb[["p.value"]][chosen_genes, 2]
actual <- hpgl_eb[["p.value"]][chosen_genes, 2]
test_that("Do cbcbSEQ and hpgltools agree on the eBayes result? (p.value, untreated)", {
  expect_equal(expected[[2]], actual[[2]], tolerance = 0.1)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 21de_limma_combat.R in ", elapsed,  " seconds.")
