start <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("066power_estimation.R:
  123\n")

small_combined <- new.env()
tt <- load(file="065_small_combined.rda", envir=small_combined)
small_combined <- small_combined[["small_combined"]]
test_proper <- sm(simple_proper(small_combined, reps=c(3,5), nsims=10))
expected <- 6
actual <- nrow(test_proper[[1]][["power_table"]])
test_that("Minimal check for proper functionality:", {
  expect_equal(expected, actual)
})

expected <- "recordedplot"
actual <- class(test_proper[[1]][["powerfd_plot"]])
test_that("Minimal check for proper plotting:", {
  expect_equal(expected, actual)
})
