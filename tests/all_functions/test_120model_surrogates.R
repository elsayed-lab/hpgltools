start <- as.POSIXlt(Sys.time())
context("120model_surrogates.R:\n")
## 2017-12, exported functions in model_surrogates:
## get_model_adjust -- note this has lots of parameters to test including:
## estimate_type: sva_supervised, sva_unsupervised, fsva, svaseq, pca,
##  ruv_supervised, ruv_residuals, ruv_empirical
## compare_surrogate_estimates,

## make_pombe_expt() invokes create_expt()
pombe_expt <- make_pombe_expt()

## Now I think on it, compare_surrogate_estimates runs all of the others for me...
pombe_result <- compare_surrogate_estimates(pombe_expt, do_catplots=TRUE)
## Hmm I am not sure what to test in this.
expected <- 4
actual <- ncol(pombe_result[["adjustments"]][["pca_adjust"]])
test_that("Do we get expected results from compare_surrogate_estimates()?", {
  expect_equal(actual, expected)
})

pombe_filt <- normalize_expt(pombe_expt, filter=TRUE)
adjust_test_sva <- get_model_adjust(pombe_filt, estimate_type="sva")
adjust_test_svaseq <- get_model_adjust(pombe_filt, estimate_type="svaseq")
adjust_test_isva <- get_model_adjust(pombe_filt, estimate_type="isva")
comparison <- cor(adjust_test_sva[["new_counts"]][, 1],
                  adjust_test_svaseq[["new_counts"]][, 1])
test_that("Do we get similar sva/svaseq results?", {
  expect_equal(comparison, 1)
})
comparison <- cor(adjust_test_sva[["new_counts"]][, 1],
                  adjust_test_isva[["new_counts"]][, 1])
test_that("Do we get similar sva/isva results?", {
  expect_gt(comparison, 0.98)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 120model_surrogates.R in ", elapsed,  " seconds."))
