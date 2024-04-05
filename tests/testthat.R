#!/usr/bin/env Rscript
library(testthat)

## test_package() is (I think) required for covr to work.
## However, when I am working on it, I want to use test_local()
## so that I can analyze the various files that are created.

test_package("hpgltools", reporter = "summary", stop_on_failure = FALSE)
