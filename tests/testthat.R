#!/usr/bin/env Rscript
<<<<<<< HEAD

library(testthat)
library(hpgltools)

print(getwd())
test_dir("tests/testthat")
=======
library(methods)
library(testthat)
library(hpgltools)

message("Beginning test_dir('testthat')")
test_result <- try(test_dir("testthat"))
if (class(test_result) == "try-error") {
    message("The test suite returned an error.")
} else {
    message("The test suite returned successfully.")
    knitr::kable(as.data.frame(test_result))
}
>>>>>>> 3d1c7f4094fa17124a141b2aeb2406119656ec68
