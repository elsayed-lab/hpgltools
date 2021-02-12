start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("125model_testing.R:
  123\n")
## 2017-12, exported functions in model_testing:
## model_test

## model_test just uses qr() to test that the data for a given statistical model
## is of sufficient rank and yells at you if it does not.
pombe_expt <- make_pombe_expt()
design <- pData(pombe_expt)

pombe_test <- model_test(design, goal = "condition")
expected <- 1
actual <- pombe_test[["condition"]]
test_that("Do we get expected rank comparisons for the pombe experiment (condition only)?", {
  expect_equal(actual, expected)
})

pombe_test <- model_test(design, goal = "batch", factors = "strain")
expected <- 1
actual <- pombe_test[["strain"]]
test_that("Do we get expected rank comparisons for the pombe experiment (~ 0 + strain + batch)?", {
  expect_equal(actual, expected)
})

pombe_test <- model_test(design, factors = "minute")
expected <- 0
actual <- pombe_test[["minute"]]
test_that("Do we get expected rank comparisons for the pombe experiment (~ 0 + strain + minute)?", {
  expect_equal(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end) - as.numeric(start))
message(paste0("\nFinished 125model_testing.R in ", elapsed,  " seconds."))
