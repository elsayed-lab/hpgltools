#!/usr/bin/env Rscript
library(methods)
library(testthat)
library(hpgltools)

all <- data.frame()
result <- 0 
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  message("Beginning test_dir('all_functions')")
  all <- try(testthat::test_dir("tests/all_functions", reporter="summary"))
  if (class(all) == "try-error") {
    result <- result + 1
    all <- data.frame()
  }
}

message("Beginning test_dir('travis')")
test_result <- try(testthat::test_dir("tests/travis", reporter="summary"))
if (class(test_result) == "try-error") {
  result <- result + 1
  test_result <- data.frame()
}

notravis <- data.frame()
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  message("Beginning test_dir('slow_tests')")
  notravis <- try(testthat::test_dir("tests/slow_tests", reporter="summary"))
  if (class(notravis) == "try-error") {
    result <- result + 1
    notravis <- data.frame()
  }
}
