#!/usr/bin/env Rscript
library(testthat)
library(hpgltools)

##test_package("hpgltools", reporter="summary", stop_on_failure=FALSE)
test_local(path = ".", reporter = "summary", stop_on_failure = FALSE)

