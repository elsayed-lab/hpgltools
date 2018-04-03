start <- as.POSIXlt(Sys.time())
context("120model_surrogates.R:\n")
## 2017-12, exported functions in model_surrogates:
## get_model_adjust -- note this has lots of parameters to test including:
## estimate_type: sva_supervised, sva_unsupervised, fsva, svaseq, pca,
##  ruv_supervised, ruv_residuals, ruv_empirical
## compare_surrogate_estimates,

## make_pombe_expt() invokes create_expt()
pombe_expt <- sm(make_pombe_expt())

## Now I think on it, compare_surrogate_estimates runs all of the others for me...
pombe_result <- sm(compare_surrogate_estimates(pombe_expt, do_catplots=TRUE))
## Hmm I am not sure what to test in this.
expected <- 4
actual <- ncol(pombe_result[["adjustments"]][["pca_adjust"]])
test_that("Do we get expected results from compare_surrogate_estimates()?", {
  expect_equal(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 120model_surrogates.R in ", elapsed,  " seconds."))
