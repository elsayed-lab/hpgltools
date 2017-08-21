#!/usr/bin/env Rscript
library(methods)
library(testthat)
library(hpgltools)

message("Beginning test_dir('testthat')")
test_result <- try(test_dir("tests/testthat"))
if (class(test_result) == "try-error") {
    message("The test suite returned an error.")
} else {
    message("The test suite returned successfully.")
    summary(as.data.frame(test_result))
}
warnings()
