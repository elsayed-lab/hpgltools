#!/usr/bin/env Rscript
library(testthat)
library(hpgltools)

test_package("hpgltools", reporter="summary")
