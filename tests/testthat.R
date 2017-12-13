#!/usr/bin/env Rscript
library(methods)
library(testthat)
library(hpgltools)

message("Beginning test_dir('travis')")
test_result <- try(testthat::test_dir("tests/travis"))
result <- 0
if (class(test_result) == "try-error") {
  result <- result + 1
  test_result <- data.frame()
}

all_functions <- data.frame()
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  message("Beginning test_dir('all_functions')")
  all <- try(testthat::test_dir("tests/all_functions"))
  if (class(all) == "try-error") {
    result <- result + 1
    all <- data.frame()
  }
}

notravis <- data.frame()
if (!identical(Sys.getenv("TRAVIS"), "true")) {
  message("Beginning test_dir('slow_tests')")
  notravis <- try(testthat::test_dir("tests/slow_tests"))
  if (class(notravis) == "try-error") {
    result <- result + 1
    notravis <- data.frame()
  }
}

if (result > 0) {
  message(paste0("There were ", result, " errors."))
}

summary(as.data.frame(test_result))
summary(as.data.frame(all))
summary(as.data.frame(notravis))
