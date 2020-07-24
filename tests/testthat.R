#!/usr/bin/env Rscript
library(testthat)
library(hpgltools)

test_check("hpgltools", reporter="summary")
