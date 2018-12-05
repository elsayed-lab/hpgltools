start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("150normalize_batch.R:\n")

pombe_expt <- make_pombe_expt()

testing <- normalize_expt(pombe_expt, filter=TRUE, batch="limma")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(27.21933, 33.08824, 147.85699, 19.53751, 89.46974,
              206.14546, 39.54991, 117.83365, 139.79517, 20.77142)
test_that("limma batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

########### THIS IS THE SOLUTION TO MY CONFUSION RECENTLY!!!!
################################                                    vvvvvvvvvvvv
testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="svaseq", surrogates=1)
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(12.7446437, 25.2166087, 119.5960057, -0.1616566, 72.0400874,
              147.9003234, 21.8864523, 84.9544863, 115.7224050, 12.3342972)
test_that("svaseq batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="fsva", surrogates=1)
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(22.62965, 36.77725, 154.51700, 18.85615, 90.74690,
              183.42573, 48.49662, 104.70992, 150.66134, 21.89822)
test_that("fsva batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

t2 <- all_adjusters(pombe_expt, filter=TRUE, estimate_type="fsva", surrogates=1)
t2$model_adjust
t2$new_counts[1:10, 1]
t2$source_counts[1:10, 1]

t3 <- counts_from_surrogates(pombe_expt, adjust=t2$model_adjust)
t3[1:10, 1]


testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="ssva")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(5.545389, 4.353055, 7.317340, 3.576433, 4.373455,
              16.194469, 6.299284, 4.905379, 7.350800, 2.123014)
test_that("ssva batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="ruvg", surrogates=1, thresh=1)
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(23, 26, 90, 32, 70, 141, 44, 64, 74, 11)
test_that("ruvg batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="ruv")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(1.93986622, 1.82724574, 7.77781751, -0.09876664, 4.99048731,
              10.05078253, 2.48853435, 5.87843369, 7.64191549, 0.80980365)
test_that("ruv batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

message("FIXME! I think something is wrong with my isva invocation.")
testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="isva")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(-6.171629, 15.531229, -13.421936, 7.192514, -20.826294,
              -78.238150, 77.874809, -88.632286, 3.597047, -11.672935)
test_that("isva batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="combat")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(20.72062, 32.99341, 150.40940, 12.93767, 87.57504,
              177.59722, 42.29628, 99.52522, 146.12594, 20.49443)
test_that("combat batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

## Now I think on it, compare_surrogate_estimates runs all of the others for me...
pombe_result <- compare_surrogate_estimates(pombe_expt, do_catplots=TRUE)
## Hmm I am not sure what to test in this.
expected <- 4
actual <- ncol(pombe_result[["adjustments"]][["pca_adjust"]])
## 01
test_that("Do we get expected results from compare_surrogate_estimates()?", {
  expect_equal(actual, expected)
})

pombe_filt <- normalize_expt(pombe_expt, filter=TRUE)
adjust_test_sva <- all_adjusters(pombe_filt, estimate_type="ssva")
adjust_test_svaseq <- all_adjusters(pombe_filt, estimate_type="svaseq")
comparison <- cor(adjust_test_sva[["new_counts"]][, 1],
                  adjust_test_svaseq[["new_counts"]][, 1])
## 02
test_that("Do we get similar sva/svaseq results?", {
  expect_gt(comparison, 0.99)
})

## Something is messed up here.
adjust_test_isva <- sm(all_adjusters(pombe_filt, estimate_type="isva"))
comparison <- cor(adjust_test_sva[["new_counts"]][, 1],
                  adjust_test_isva[["new_counts"]][, 1], method="spearman")
## 03
test_that("Do we get similar sva/isva results?", {
  expect_gt(comparison, 0.55)
})

adjust_test_smartsva <- all_adjusters(pombe_filt, estimate_type="smartsva")
comparison <- cor(adjust_test_smartsva[["new_counts"]][, 1],
                  adjust_test_sva[["new_counts"]][, 1])
## 04
test_that("Do we get similar sva/isva results?", {
  expect_gt(comparison, 0.99)
})

adjust_test_ruv <- all_adjusters(pombe_filt, estimate_type="ruv")
comparison <- cor(adjust_test_ruv[["new_counts"]][, 1],
                  adjust_test_sva[["new_counts"]][, 1])
## 05
test_that("Do we get similar sva/isva results?", {
  expect_gt(comparison, 0.99)
})

## If I invoke RUVg, it returns surrogates + counts.
## I attempted to implement their solver, I therefore want to test that
## a matrix acquired from their solver + their model gets similar/identical
## counts to that acquired entirely from RUVg.  If this is true, then I think
## that gives me some confidence that I can use the solver for other methods.
ruv_estimate <- all_adjusters(pombe_filt, estimate_type="ruvg")
ruv_counts <- ruv_estimate$source_counts
ruv_model <- ruv_estimate$model_adjust
my_counts <- counts_from_surrogates(pombe_filt, adjust=ruv_model, method="ruv")
test_that("Do we get to a similar end point from ruv with/without solving?", {
  expect_equal(ruv_counts, my_counts, tolerance=0.5)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 120model_surrogates.R in ", elapsed,  " seconds."))
