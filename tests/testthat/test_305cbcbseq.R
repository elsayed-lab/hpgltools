start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
library(pasilla)
tt <- sm(library(edgeR))
## data(pasillaGenes)
num_installed <- please_install("kokrah/cbcbSEQ")
cbcb <- sm(library(cbcbSEQ))
context("05cbcbseq.R: Compare cbcbSEQ output to hpgltools.\n")
## This test is intended to compare Kwame Okrah/Hector Corrada Bravo's cbcbSEQ output
## against that received from this.  hpgltools is a derivative of it, therefore it should
## provide identical results when called in a similar fashion.

## These are being performed as per the cbcbSEQ vignette:
## 1. Reads were read manually from the pasilla dataset
## 2. log2(quantile(cpm(data))) was performed
## 3. A slightly modified version of voom() was performed
## 4. lmfit()->topTables() was performed
## Hopefully hpgltools can recapitulate the results.

### TODO:
## I recently made explicit the reading of annotation and read data into the expressionset.
## As a result, the order of rows in the count tables may now be changed to match that of
## the count tables.  This results in potential gene-order changes between cbcbSEQ and hpgltools.
## Therefore many tests may fail not because the data is wrong, but because of gene-order changes.

### Notes:
## cbcbSEQ::log2CPM() returns the normalized libsizes rather than those following log2() transformation.
## It then goes on to call voom() with this libsize rather than the log2().
## In order to follow this behaivor, my normalization functions keep a copy of the libsizes
## and count tables from beginning to end of all analyses to ensure that this is still accessible.
##
## I copied the modified voom() from cbcbSEQ and made a couple of small changes to make it more flexible
## vis a vis data which has been log transformed and/or converted via cpm()
## This may be a sticking point

## Read the raw pasilla data into a dataframe and make an expressionset
datafile <- system.file("extdata/pasilla_gene_counts.tsv", package = "pasilla")
## Load the counts and drop super-low counts genes
counts <- read.table(datafile, header = TRUE, row.names = 1)
counts <- counts[rowSums(counts) > ncol(counts), ]
## Set up a quick design to be used by cbcbSEQ and hpgltools
design <- data.frame(row.names = colnames(counts),
                     condition = c("untreated","untreated","untreated",
                                 "untreated","treated","treated","treated"),
                     libType = c("single_end","single_end","paired_end",
                               "paired_end","single_end","paired_end","paired_end"))
metadata <- design
colnames(metadata) <- c("condition", "batch")
metadata[["sampleid"]] <- rownames(metadata)
cbcb_data <- as.matrix(counts)

pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]
hpgl_data <- exprs(pasilla_expt)

## Check that normalization tools work similarly
cbcb_quantile <- cbcbSEQ::qNorm(cbcb_data)
hpgl_quantile_data <- normalize_expt(pasilla_expt, transform = "raw", norm = "quant",
                                     convert = "raw", filter = FALSE)
hpgl_quantile <- exprs(hpgl_quantile_data)
expected <- cbcb_quantile[sort(rownames(cbcb_quantile)), ]
actual <- hpgl_quantile[sort(rownames(hpgl_quantile)), ]

test_that("Are the quantile normalizations identical?", {
  expect_equal(expected, actual)
})

## Check that quantile() normalizations are identical
cbcb_qcpm <- cbcbSEQ::qNorm(cbcb_data)
cbcb_edger_qcpm <- edgeR::cpm(cbcb_qcpm)
cbcb_quantile <- cbcbSEQ::qNorm(cbcb_data)
expected <- cbcb_quantile[sort(rownames(cbcb_quantile)), ]
hpgl_quantile <- hpgl_norm(pasilla_expt, norm = "quant")
hpgl_quantile <- hpgl_quantile[["count_table"]]
actual <- hpgl_quantile[sort(rownames(hpgl_quantile)), ]
test_that("Are quantile() normalizations identical?", {
  expect_equal(expected, actual)
})

## Check that cpm(quantile()) normalizations are identical
hpgl_qcpm <- hpgl_norm(pasilla_expt, norm = "quant", convert = "cpm", filter = FALSE)
hpgl_qcpm <- hpgl_qcpm[["count_table"]]
hpgl_qcpm <- hpgl_qcpm[sort(rownames(hpgl_qcpm)), ]
expected <- cbcb_edger_qcpm
actual <- hpgl_qcpm
test_that("Are cpm(quantile()) conversions identical?", {
  expect_equal(expected, actual)
})

## Check that log2(cpm(quantile())) are identical
## There are a couple of ways to invoke these normalizations, so I will
## explicitly (and hopefully redundantly) invoke both
cbcb_l2qcpm_data <- cbcbSEQ::log2CPM(cbcb_quantile)
cbcb_l2qcpm <- cbcb_l2qcpm_data[["y"]]
hpgl_l2qcpm_data <- hpgl_norm(pasilla_expt, transform = "log2", norm = "quant",
                              convert = "cbcbcpm", filter = FALSE)
hpgl_l2qcpm <- hpgl_l2qcpm_data[["count_table"]]
hpgl_l2qcpm <- hpgl_l2qcpm[sort(rownames(hpgl_l2qcpm)), ]
hpgl_l2qcpm_expt <- normalize_expt(pasilla_expt, transform = "log2", norm = "quant",
                                   convert = "cbcbcpm", filter = FALSE)
hpgl_l2qcpm2 <- exprs(hpgl_l2qcpm_expt)
hpgl_l2qcpm2 <- hpgl_l2qcpm2[sort(rownames(hpgl_l2qcpm2)), ]
expected <- cbcb_l2qcpm
actual <- hpgl_l2qcpm
test_that("Are l2qcpm conversions/transformations identical using cbcbSEQ vs. hpgl_norm()?", {
  expect_equal(expected, actual)
})
actual <- hpgl_l2qcpm2
test_that("Are l2qcpm conversions/transformations identical using cbcbSEQ vs. normalize_expt()?", {
  expect_equal(expected, actual)
})

## Check that the libsizes are properly maintained
cbcb_libsize <- cbcb_l2qcpm_data[["lib.size"]]
hpgl_libsize <- hpgl_l2qcpm_data[["intermediate_counts"]][["normalization"]][["libsize"]]
expected <- cbcb_libsize
actual <- hpgl_libsize
test_that("In preparing for voom(), are the library sizes maintained?", {
  expect_equal(expected, actual)
})

## If we get here without problems, then voom->topTable should be ready to go without problems.
## Given that, lets try a voom() invocation and see what happens.
condition <- design[["condition"]]
test_model <- model.matrix(~condition)
cbcb_voom <- suppressWarnings(cbcbSEQ::voomMod(x = as.matrix(cbcb_l2qcpm), design = test_model, lib.size = cbcb_libsize))
hpgl_voom <- suppressWarnings(cbcbSEQ::voomMod(x = as.matrix(hpgl_l2qcpm), design = test_model, lib.size = hpgl_libsize))
hpgl_voom2 <- suppressWarnings(hpgltools::hpgl_voom(as.matrix(hpgl_l2qcpm), model = test_model,
                                   libsize = hpgl_libsize, logged = TRUE, converted = TRUE))
hpgl_voom3 <- suppressWarnings(hpgltools::hpgl_voom(as.matrix(hpgl_quantile), test_model,
                                                    libsize = hpgl_libsize, logged = FALSE, converted = FALSE))
expected <- cbcb_voom
actual <- hpgl_voom
test_that("Do different voom() invocations end with the same result?", {
  expect_equal(expected, actual)
})

expected <- cbcb_voom[["E"]]
actual <- hpgl_voom2[["E"]]
test_that("Does calling cbcbSEQ voom with hpgl-modified data return the same result?", {
  expect_equal(expected, actual)
})
actual <- hpgl_voom3[["E"]]
test_that("Does calling hpgltools::voom with hpgl-modified data return the same result?", {
  expect_equal(expected, actual)
})

## To be extra-paranoid, make sure that the limma_pairwise() function invokes voom correctly.
## Note that this is where the data-ordering problems appear.
hpgl_limma <- suppressWarnings(limma_pairwise(input = hpgl_l2qcpm_expt, model_batch = FALSE, limma_method = "ls",
                                              model_intercept = TRUE, which_voom = "hpgl"))

## First check the voom result from limma_pairwise
hpgl_limma_voom <- hpgl_limma[["voom_result"]]
hpgl_limma_voom_e <- hpgl_limma[["voom_result"]][["E"]][order(rownames(hpgl_limma[["voom_result"]][["E"]])), ]
cbcb_voom_e <- cbcb_voom[["E"]][order(rownames(cbcb_voom[["E"]])), ]
test_that("Limma results, voom.", {
  expect_equal(cbcb_voom_e, hpgl_limma_voom_e)
})

## Then the result from lmFit
hpgl_limma_fit_coef <- hpgl_limma[["fit"]][["coefficients"]][order(rownames(hpgl_limma[["fit"]][["coefficients"]])), ]
cbcb_fit <- limma::lmFit(cbcb_voom)
cbcb_fit_coef <- cbcb_fit[["coefficients"]][order(rownames(cbcb_fit[["coefficients"]])), ]
colnames(cbcb_fit_coef) <- c("(Intercept)", "untreated")
test_that("Limma results, fitting coefficients.", {
  expect_equal(cbcb_fit_coef, hpgl_limma_fit_coef)
})

hpgl_limma_fit_stdev <- hpgl_limma[["fit"]][["stdev.unscaled"]][order(rownames(hpgl_limma[["fit"]][["stdev.unscaled"]])), ]
cbcb_fit_stdev <- cbcb_fit[["stdev.unscaled"]][order(rownames(cbcb_fit[["stdev.unscaled"]])), ]
colnames(cbcb_fit_stdev) <- c("(Intercept)", "untreated")
test_that("Limma results, fitting standard deviations.", {
  expect_equal(cbcb_fit_stdev, hpgl_limma_fit_stdev)
})

## Now the result from eBayes
cbcb_eb <- limma::eBayes(cbcb_fit)
hpgl_eb <- limma::eBayes(hpgl_limma[["fit"]])
test_that("Limma results, eBayes.", {
  expect_equal(sort(cbcb_eb[["F"]]), sort(hpgl_eb[["F"]]))
})

cbcb_top <- limma::topTable(cbcb_eb, number = nrow(cbcb_eb))
cbcb_top <- cbcb_top[sort(rownames(cbcb_top)), ]
hpgl_top <- limma::topTable(hpgl_eb, number = nrow(cbcb_eb))
hpgl_top <- hpgl_top[sort(rownames(hpgl_top)), ]
test_that("Limma results, toptable.", {
  expect_equal(cbcb_top, hpgl_top)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 05cbcbseq.R in ", elapsed, " seconds.")
