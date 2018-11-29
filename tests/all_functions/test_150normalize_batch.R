start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("150normalize_batch.R:\n")
## 2017-12, exported functions in annotation_gff:

pombe_expt <- make_pombe_expt()

testing <- normalize_expt(pombe_expt, filter=TRUE, batch="limma")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(27.21933, 33.08824, 147.85699, 19.53751, 89.46974,
              206.14546, 39.54991, 117.83365, 139.79517, 20.77142)
test_that("limma batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="svaseq")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(12.7446437, 25.2166087, 119.5960057, -0.1616566, 72.0400874,
              147.9003234, 21.8864523, 84.9544863, 115.7224050, 12.3342972)
test_that("svaseq batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="fsva")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(22.62965, 36.77725, 154.51700, 18.85615, 90.74690,
              183.42573, 48.49662, 104.70992, 150.66134, 21.89822)
test_that("fsva batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="ssva")
test_counts <- exprs(testing)
actual <- as.numeric(test_counts[1:10, 1])
expected <- c(5.545389, 4.353055, 7.317340, 3.576433, 4.373455,
              16.194469, 6.299284, 4.905379, 7.350800, 2.123014)
test_that("ssva batch modification provides expected values?", {
  expect_equal(expected, actual, tolerance=0.0001)
})

testing <-  normalize_expt(pombe_expt, filter=TRUE, batch="ruvg")
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

end <- as.POSIXlt(Sys.time())
elapsed <- round(x=as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 150normalize_batch.R in ", elapsed,  " seconds."))
