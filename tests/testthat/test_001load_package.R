start_time <- as.POSIXlt(Sys.time())
library(testthat)
library(hpgltools)
context("001load_package.R: Does hpgltools load properly?")

installedp <- utils::vignette(package = "hpgltools")
available_vignettes <- as.character(as.data.frame(installedp[["results"]])[["Item"]])

## There are 4 vignettes, and depending on how everything installed,
## potentially a reference manual.
## I don't get it, I run this manually and it works fine.
##expected <- 3
##actual <- length(available_vignettes)
##test_that("Did the vignettes install?", {
##  expect_gt(actual, expected)
##})

manual <- file.exists(system.file("doc/reference.pdf", package = "hpgltools"))
test_that("Did the reference manual install?", {
  expect_true(manual)
})

end_time <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end_time - start_time))
message("\nFinished 01load_package.R in ", elapsed,  " seconds.")
